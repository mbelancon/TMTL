"""Servidor local (somente biblioteca padrão) que expõe o motor de ajuste à interface web."""
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
from .models import MODELS, Accuracy
from .optim import FitProblem, Stopped, make_free, run_fit
from .batch import BatchJob

STATIC = Path(__file__).parent / "static"
JOBS: dict[str, "Job"] = {}
# vários usuários no mesmo servidor: limita ajustes simultâneos (os demais esperam na fila)
FIT_SLOTS = threading.BoundedSemaphore(int(__import__("os").environ.get("TLFIT_MAX_FITS", "2")))
JOB_TTL = 3600  # s: ajustes terminados são descartados depois de 1 h
LAST_PING = [time.time()]  # modo aplicativo: encerra quando a janela fecha


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
    def __init__(self, prob: FitProblem, opts: dict):
        self.id = uuid.uuid4().hex[:10]
        self.prob = prob
        self.opts = opts
        self.stop = threading.Event()
        self.lock = threading.Lock()
        self.status = "running"
        self.phase = "pso"
        self.iteration = 0
        self.run_idx, self.n_runs = 0, 1
        self.n_iter = int(opts.get("iterations", 40))
        self.history = []          # [n_eval, melhor χ², fase]
        self.swarm = None          # posições atuais (unidades físicas)
        self.swarm_chi = None
        self.trail = []            # posições das últimas iterações
        self.log = []
        self.result = None
        self.started = time.time()
        self._last_hist = 0.0

    def say(self, msg):
        with self.lock:
            self.log.append(f"[{time.time() - self.started:6.1f} s] {msg}")

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
            if info["phase"] == "grade":
                # varredura em dois níveis: mostra todos os pontos já calculados
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
                    self.trail = []          # nova rodada: limpa o rastro do enxame
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
                self.history.append([self.prob.n_eval, res["chi2"], "fim"])
        except Stopped:
            self.say("Ajuste interrompido pelo usuário; mantendo o melhor ponto encontrado.")
            with self.lock:
                self.status = "stopped"
                if self.prob.best_z is not None:
                    z = self.prob.best_z
                    self.result = dict(params=self.prob.phys(z), chi2=self.prob.best_chi,
                                       sigma={}, extra=None, curve=self.prob.best_curve,
                                       n_eval=self.prob.n_eval,
                                       seconds=time.time() - self.started)
        except Exception as e:  # noqa: BLE001
            self.say("ERRO: " + "".join(traceback.format_exception_only(type(e), e)).strip())
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
    out = []
    for m in MODELS.values():
        out.append(dict(key=m.key, title=m.title, description=m.description, presets=m.presets,
                        linear=list(m.linear), theta_info=m.theta_info,
                        proportional=sorted(set(m.linear) | set(getattr(m, "proportional", ()))),
                        params=[dict(name=p.name, label=p.label, unit=p.unit, value=p.value,
                                     fit=p.fit, lo=p.lo, hi=p.hi) for p in m.params]))
    return out


def _acc(level):
    return Accuracy().scaled(float(level or 1.0))


class Handler(BaseHTTPRequestHandler):
    server_version = "tlfit/1.0"

    def log_message(self, fmt, *args):  # silencioso
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

    def do_GET(self):
        path = self.path.split("?")[0]
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
                # linha opcional "# tlfit: modelo=<chave> m=… w0=…": escolhe o modelo e a geometria
                hdr = re.search(r"^#\s*tlfit:\s*modelo=(\w+)(.*)$", text, re.M)
                extra = {}
                if hdr and hdr.group(1) in MODELS:
                    extra = dict(model=hdr.group(1), params={
                        k: float(v) for k, v in re.findall(r"(\w+)=([-+0-9.eE]+)(?:\s|$)", hdr.group(2))})
                    pr = re.search(r"preset=(\S+)", hdr.group(2))
                    if pr:   # os arquivos de exemplo são ASCII: "Agua" casa com "Água"
                        want = _ascii(pr.group(1)).lower()
                        extra["preset"] = next((k for k in MODELS[hdr.group(1)].presets
                                                if _ascii(k).lower() == want), pr.group(1))
                return self._send(200, dict(name=name, text=text, **extra))
            return self._send(404, {"error": "exemplo não encontrado"})
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
                return self._send(404, {"error": "série desconhecida"})
            return self._send(200, job.snapshot())
        if path.startswith("/api/job/"):
            job = JOBS.get(path.rsplit("/", 1)[-1])
            if not job:
                return self._send(404, {"error": "job desconhecido"})
            return self._send(200, job.snapshot())
        self._send(404, {"error": "não encontrado"})

    def do_POST(self):
        path = self.path.split("?")[0]
        try:
            req = self._json()
            if path == "/api/parse":
                data, header = parse_table(req["text"])
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
                    return self._send(400, {"error": "marque ao menos um parâmetro para ajustar"})
                for f in free:  # valor inicial dentro dos limites
                    p[f.name] = float(np.clip(p[f.name], f.lo, f.hi))
                prob = FitProblem(m, req["t"], req["y"], p, free, _acc(req.get("accuracy")))
                _prune()
                job = Job(prob, req.get("options", {}))
                JOBS[job.id] = job
                threading.Thread(target=job.run, daemon=True).start()
                return self._send(200, {"id": job.id})
            if path == "/api/batch":
                m = MODELS[req["model"]]
                p = m.defaults()
                p.update({k: float(v) for k, v in req["params"].items()})
                if not req["free"]:
                    return self._send(400, {"error": "marque ao menos um parâmetro para ajustar"})
                items = [it for it in req["items"] if len(it.get("t", [])) >= 3]
                if not items:
                    return self._send(400, {"error": "nenhum transiente com pelo menos 3 pontos"})
                _prune()
                job = BatchJob(m, p, req["free"], items, req.get("options", {}),
                               _acc(req.get("accuracy")), slots=FIT_SLOTS)
                JOBS[job.id] = job
                threading.Thread(target=job.run, daemon=True).start()
                return self._send(200, {"id": job.id})
            if (path.startswith("/api/job/") or path.startswith("/api/batch/")) and path.endswith("/stop"):
                job = JOBS.get(path.split("/")[3])
                if job:
                    job.stop.set()
                return self._send(200, {"ok": True})
            self._send(404, {"error": "não encontrado"})
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            self._send(400, {"error": str(e)})


def make_server(port=0, host="127.0.0.1"):
    """Cria o servidor sem iniciá-lo (usado pelo aplicativo de desktop)."""
    return ThreadingHTTPServer((host, port), Handler)


def serve(port=8765, open_browser=True, host="127.0.0.1"):
    httpd = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{'127.0.0.1' if host in ('0.0.0.0', '') else host}:{httpd.server_address[1]}/"
    print(f"Interface de ajuste em {url}  (Ctrl+C para encerrar)")
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
