"""Local server (standard library only) that exposes the fitting engine to the web interface.

Messages are produced in English and translated to the language sent by the
interface ("lang": "en" or "pt"; see tlfit/i18n.py)."""
from __future__ import annotations

import json
import re
import unicodedata
import threading
import time
import traceback
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

from .dataio import guess_power, parse_table
from .i18n import Msg, both, norm_lang, tr
from .models import MODELS, PRESET_ALIASES, Accuracy
from .optim import FitProblem, Stopped, make_free, run_fit
from .batch import BatchJob

STATIC = Path(__file__).parent / "static"
JOBS: dict[str, "Job"] = {}
# several users on the same server: limits simultaneous fits (the others wait in the queue)
FIT_SLOTS = threading.BoundedSemaphore(int(__import__("os").environ.get("TLFIT_MAX_FITS", "2")))
JOB_TTL = 3600  # s: finished fits are discarded after 1 h
LAST_PING = [time.time()]  # application mode: exits when the window closes


def _prune():
    now = time.time()
    for k in [k for k, j in JOBS.items() if j.status != "running" and j.status != "queued"
              and now - j.started > JOB_TTL]:
        JOBS.pop(k, None)


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


class Job:
    def __init__(self, prob: FitProblem, opts: dict, lang="en"):
        self.id = uuid.uuid4().hex[:10]
        self.lang = lang
        self.prob = prob
        self.opts = opts
        self.stop = threading.Event()
        self.lock = threading.Lock()
        self.status = "running"
        self.phase = "pso"
        self.iteration = 0
        self.run_idx, self.n_runs = 0, 1
        self.n_iter = int(opts.get("iterations", 40))
        self.history = []          # [n_eval, best χ², phase]
        self.swarm = None          # current positions (physical units)
        self.swarm_chi = None
        self.trail = []            # positions of the last iterations
        self.log = []
        self.result = None
        self.started = time.time()
        self._last_hist = 0.0

    def say(self, msg):
        with self.lock:
            self.log.append(f"[{time.time() - self.started:6.1f} s] {tr(msg, self.lang)}")

    def progress(self, info):
        prob = self.prob
        with self.lock:
            self.phase = info["phase"]
            self.run_idx, self.n_runs = info.get("run", 0), info.get("n_runs", 1)
            chi = prob.best_chi
            if np.isfinite(chi):
                if not self.history or chi < self.history[-1][1] or \
                        time.time() - self._last_hist > 0.5:
                    self.history.append([prob.n_eval, chi, info["phase"]])
                    self._last_hist = time.time()
            if info["phase"] == "grid":
                # two-level sweep: shows all points already computed
                self.iteration, self.n_iter = info["iteration"], info["n_iter"]
                self.swarm = np.array([[fp.to_phys(z) for fp, z in zip(prob.free, row)]
                                       for row in info["positions"]])
                self.swarm_chi = info["chi"]
                self.trail = []
            elif info["phase"] == "pso":
                self.iteration = info["iteration"]
                pos = np.array([[fp.to_phys(z) for fp, z in zip(prob.free, row)]
                                for row in info["positions"]])
                self.swarm = pos
                self.swarm_chi = info["chi"]
                if info.get("iteration", 0) == 0:
                    self.trail = []          # new run: clears the swarm trail
                self.trail.append(pos)
                self.trail = self.trail[-12:]

    def run(self):
        self.status = "queued"
        with FIT_SLOTS:
            self.status = "running"
            self.started = time.time()
            self._run()

    def _run(self):
        try:
            res = run_fit(self.prob, method=self.opts.get("method", "auto"),
                          n_particles=int(self.opts.get("particles", 20)),
                          n_iter=self.n_iter, seed=self.opts.get("seed"),
                          restarts=int(self.opts.get("restarts", 3)),
                          callback=self.progress, stop=self.stop, log=self.say)
            with self.lock:
                self.result = res
                self.status = "done"
                self.history.append([self.prob.n_eval, res["chi2"], "end"])
        except Stopped:
            self.say("Fit stopped by the user; keeping the best point found.")
            with self.lock:
                self.status = "stopped"
                if self.prob.best_z is not None:
                    z = self.prob.best_z
                    self.result = dict(params=self.prob.phys(z), chi2=self.prob.best_chi,
                                       sigma={}, extra=None, curve=self.prob.best_curve,
                                       n_eval=self.prob.n_eval,
                                       seconds=time.time() - self.started)
        except Exception as e:  # noqa: BLE001
            self.say(Msg("ERROR: {err}", err="".join(traceback.format_exception_only(type(e), e)).strip()))
            traceback.print_exc()
            with self.lock:
                self.status = "error"

    def snapshot(self):
        prob = self.prob
        with self.lock:
            best = prob.phys(prob.best_z) if prob.best_z is not None else None
            return _clean(dict(
                id=self.id, status=self.status, phase=self.phase, iteration=self.iteration,
                run=self.run_idx, n_runs=self.n_runs,
                n_iter=self.n_iter, n_eval=prob.n_eval, best_chi=prob.best_chi,
                best=best, best_curve=prob.best_curve, history=self.history,
                free=[dict(name=f.name, lo=f.lo, hi=f.hi, log=f.log) for f in prob.free],
                swarm=self.swarm, swarm_chi=self.swarm_chi, trail=self.trail,
                log=self.log, result=self.result,
                elapsed=time.time() - self.started))


def _ascii(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()


def _model_info():
    """Models for the interface; every text goes as {"en": ..., "pt": ...}."""
    out = []
    for m in MODELS.values():
        theta = {n: {k: both(v) for k, v in d.items()} for n, d in m.theta_info.items()}
        out.append(dict(key=m.key, title=both(m.title), description=both(m.description),
                        presets=m.presets, preset_labels={k: both(k) for k in m.presets},
                        linear=list(m.linear), theta_info=theta,
                        proportional=sorted(set(m.linear) | set(getattr(m, "proportional", ()))),
                        params=[dict(name=p.name, label=both(p.label), unit=both(p.unit), value=p.value,
                                     fit=p.fit, lo=p.lo, hi=p.hi) for p in m.params]))
    return out


def _acc(level):
    return Accuracy().scaled(float(level or 1.0))


class Handler(BaseHTTPRequestHandler):
    server_version = "tlfit/1.0"

    def log_message(self, fmt, *args):  # silent
        pass

    def _send(self, code, body, ctype="application/json"):
        if not isinstance(body, (bytes, bytearray)):
            body = json.dumps(_clean(body)).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self):
        n = int(self.headers.get("Content-Length", 0))
        return json.loads(self.rfile.read(n) or b"{}")

    def _lang_q(self):
        q = self.path.split("?", 1)[1] if "?" in self.path else ""
        m = re.search(r"(?:^|&)lang=([A-Za-z_-]+)", q)
        return norm_lang(m.group(1) if m else "en")

    def do_GET(self):
        path = self.path.split("?")[0]
        lang = self._lang_q()
        if path in ("/", "/index.html"):
            return self._send(200, (STATIC / "index.html").read_bytes(), "text/html; charset=utf-8")
        if path == "/api/examples":
            ex = sorted(p.name for p in (Path(__file__).parent / "examples").glob("*.dat"))
            return self._send(200, ex)
        if path.startswith("/api/example/"):
            name = Path(path.rsplit("/", 1)[-1]).name
            f = Path(__file__).parent / "examples" / name
            if f.is_file():
                text = f.read_text(encoding="utf-8")
                # optional line "# tlfit: model=<key> m=… w0=…" (older files: "modelo="):
                # chooses the model and the geometry
                hdr = re.search(r"^#\s*tlfit:\s*model[o]?=(\w+)(.*)$", text, re.M)
                extra = {}
                if hdr and hdr.group(1) in MODELS:
                    extra = dict(model=hdr.group(1), params={
                        k: float(v) for k, v in re.findall(r"(\w+)=([-+0-9.eE]+)(?:\s|$)", hdr.group(2))})
                    pr = re.search(r"preset=(\S+)", hdr.group(2))
                    if pr:   # example files are ASCII; Portuguese names ("Agua") are aliases
                        want = _ascii(pr.group(1)).lower()
                        want = _ascii(PRESET_ALIASES.get(want, want)).lower()
                        extra["preset"] = next((k for k in MODELS[hdr.group(1)].presets
                                                if _ascii(k).lower() == want), pr.group(1))
                return self._send(200, dict(name=name, text=text, **extra))
            return self._send(404, {"error": tr("example not found", lang)})
        if path == "/api/ping":
            LAST_PING[0] = time.time()
            return self._send(200, {"ok": True})
        if path == "/api/models":
            return self._send(200, _model_info())
        if path == "/api/series-examples":
            base = Path(__file__).parent / "examples"
            out = []
            for d in sorted(p for p in base.iterdir() if p.is_dir() and p.name.startswith("serie_")):
                out.append(dict(name=d.name, files=[dict(name=f.name, text=f.read_text(encoding="utf-8", errors="replace"))
                                                    for f in sorted(d.glob("*.dat"))]))
            return self._send(200, out)
        if path.startswith("/api/batch/"):
            job = JOBS.get(path.rsplit("/", 1)[-1])
            if not job:
                return self._send(404, {"error": tr("unknown series", lang)})
            return self._send(200, job.snapshot())
        if path.startswith("/api/job/"):
            job = JOBS.get(path.rsplit("/", 1)[-1])
            if not job:
                return self._send(404, {"error": tr("unknown job", lang)})
            return self._send(200, job.snapshot())
        self._send(404, {"error": tr("not found", lang)})

    def do_POST(self):
        path = self.path.split("?")[0]
        lang = self._lang_q()
        try:
            req = self._json()
            lang = norm_lang(req.get("lang", lang))
            if path == "/api/parse":
                data, header = parse_table(req["text"])
                header = [tr(h, lang) for h in header]
                power, src = guess_power(req.get("name", ""), req["text"])
                return self._send(200, dict(data=data.T.tolist(), header=header,
                                            power=power, power_source=src))
            if path == "/api/simulate":
                m = MODELS[req["model"]]
                p = m.defaults()
                p.update({k: float(v) for k, v in req["params"].items()})
                t = np.asarray(req["t"], float)
                t0 = time.time()
                info = {}
                y = m.intensity(t, p, _acc(req.get("accuracy")), info=info)
                return self._send(200, dict(y=y, ms=1000 * (time.time() - t0), info=info))
            if path == "/api/fit":
                m = MODELS[req["model"]]
                p = m.defaults()
                p.update({k: float(v) for k, v in req["params"].items()})
                free = [make_free(f["name"], f["lo"], f["hi"], f.get("log")) for f in req["free"]]
                if not free:
                    return self._send(400, {"error": tr("select at least one parameter to fit", lang)})
                for f in free:  # initial value within the bounds
                    p[f.name] = float(np.clip(p[f.name], f.lo, f.hi))
                prob = FitProblem(m, req["t"], req["y"], p, free, _acc(req.get("accuracy")))
                _prune()
                job = Job(prob, req.get("options", {}), lang)
                JOBS[job.id] = job
                threading.Thread(target=job.run, daemon=True).start()
                return self._send(200, {"id": job.id})
            if path == "/api/batch":
                m = MODELS[req["model"]]
                p = m.defaults()
                p.update({k: float(v) for k, v in req["params"].items()})
                if not req["free"]:
                    return self._send(400, {"error": tr("select at least one parameter to fit", lang)})
                items = [it for it in req["items"] if len(it.get("t", [])) >= 3]
                if not items:
                    return self._send(400, {"error": tr("no transient with at least 3 points", lang)})
                _prune()
                job = BatchJob(m, p, req["free"], items, req.get("options", {}),
                               _acc(req.get("accuracy")), slots=FIT_SLOTS, lang=lang)
                JOBS[job.id] = job
                threading.Thread(target=job.run, daemon=True).start()
                return self._send(200, {"id": job.id})
            if (path.startswith("/api/job/") or path.startswith("/api/batch/")) and path.endswith("/stop"):
                job = JOBS.get(path.split("/")[3])
                if job:
                    job.stop.set()
                return self._send(200, {"ok": True})
            self._send(404, {"error": tr("not found", lang)})
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            self._send(400, {"error": tr(e.args[0], lang) if e.args else str(e)})


def make_server(port=0, host="127.0.0.1"):
    """Creates the server without starting it (used by the desktop application)."""
    return ThreadingHTTPServer((host, port), Handler)


def serve(port=8765, open_browser=True, host="127.0.0.1"):
    httpd = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{'127.0.0.1' if host in ('0.0.0.0', '') else host}:{httpd.server_address[1]}/"
    print(f"Fitting interface at {url}  (Ctrl+C to quit)")
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
