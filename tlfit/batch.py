"""Ajuste de uma série de transientes (por exemplo, várias potências de excitação)
e a estatística dos parâmetros ao longo da série.

Análises usuais em lente/espelho térmico:
* θ ∝ P_e: θ = −P_e A_e L (ds/dT) φ /(K λ_p). Ajusta-se θ × P por uma reta pela
  origem (inclinação θ/P_e) e também com intercepto, para checar a linearidade.
* D (ou t_c) não deve depender de P: média, desvio-padrão e a inclinação de D × P,
  que acusa aquecimento ou outros efeitos quando é significativa.
"""
from __future__ import annotations

import threading
import time
import traceback
import uuid

import numpy as np

from .models import Model
from .optim import FitProblem, Stopped, make_free, run_fit


def fit_series(model: Model, items, params: dict, free_specs, options=None, acc=None,
               callback=None, stop=None, log=print, on_result=None):
    """items: lista de dicts {name, t, y, x}. Ajusta cada um com as mesmas
    configurações e devolve a lista de resultados (na ordem de entrada)."""
    options = dict(options or {})
    chain = options.get("chain", False)
    method = options.get("method", "auto")
    results = []
    prev = None
    n = len(items)
    for k, it in enumerate(items):
        if stop is not None and stop.is_set():
            raise Stopped
        t0 = time.time()
        p = dict(params)
        if chain and prev is not None:
            p.update({f["name"]: prev[f["name"]] for f in free_specs})
        free = [make_free(f["name"], f["lo"], f["hi"], f.get("log")) for f in free_specs]
        for f in free:
            p[f.name] = float(np.clip(p[f.name], f.lo, f.hi))
        kw = {} if acc is None else {"acc": acc}
        prob = FitProblem(model, it["t"], it["y"], p, free, **kw)

        def cb(info, k=k):
            if callback:
                callback(dict(info, item=k, n_items=n))

        try:
            r = run_fit(prob, method=method, n_particles=int(options.get("particles", 20)),
                        n_iter=int(options.get("iterations", 30)),
                        restarts=int(options.get("restarts", 3)), callback=cb, stop=stop,
                        log=lambda s: None)
        except Stopped:
            raise
        except Exception as e:  # noqa: BLE001 — um arquivo ruim não derruba a série
            log(f"{it['name']}: falhou ({e})")
            res = dict(name=it["name"], x=it.get("x"), status="erro", error=str(e))
            results.append(res)
            if on_result:
                on_result(k, res)
            continue
        y = np.asarray(it["y"], float)
        ss_res = float(((r["curve"] - y) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        extra = r.get("extra") or {}
        res = dict(name=it["name"], x=it.get("x"), status="ok",
                   params={f["name"]: r["params"][f["name"]] for f in free_specs},
                   sigma=r.get("sigma") or {}, chi2=r["chi2"],
                   r2=1 - ss_res / ss_tot if ss_tot > 0 else None,
                   at_bound=extra.get("at_bound", []), t=list(map(float, it["t"])),
                   curve=r["curve"].tolist(), seconds=time.time() - t0)
        results.append(res)
        if on_result:
            on_result(k, res)
        prev = r["params"]
        vals = ", ".join(f"{nm} = {v:.5g}" for nm, v in res["params"].items())
        warn = f"  (no limite: {', '.join(res['at_bound'])})" if res["at_bound"] else ""
        log(f"{k + 1}/{n} {it['name']}: {vals}; χ² = {r['chi2']:.3e} "
            f"({res['seconds']:.1f} s){warn}")
    return results


def _line_origin(x, y):
    """y = a·x. Erro de a pela dispersão dos pontos em torno da reta."""
    sxx = float(x @ x)
    a = float(x @ y) / sxx
    n = len(x)
    s2 = float(((y - a * x) ** 2).sum()) / max(1, n - 1)
    return a, float(np.sqrt(s2 / sxx)) if n > 1 else float("nan")


def _line(x, y):
    """y = a·x + b, mínimos quadrados ordinários, com erros e R²."""
    n = len(x)
    if n < 3:
        return None
    A = np.c_[x, np.ones(n)]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    s2 = float(res @ res) / (n - 2)
    cov = s2 * np.linalg.inv(A.T @ A)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return dict(slope=float(coef[0]), slope_err=float(np.sqrt(cov[0, 0])),
                intercept=float(coef[1]), intercept_err=float(np.sqrt(cov[1, 1])),
                r2=1 - float(res @ res) / ss_tot if ss_tot > 0 else None)


def series_stats(model: Model, results):
    """Estatística de cada parâmetro ajustado ao longo da série."""
    ok = [r for r in results if r.get("status") == "ok" and r.get("x") is not None]
    if not ok:
        return {}
    names = list(ok[0]["params"].keys())
    x = np.array([r["x"] for r in ok], float)
    out = {}
    for nm in names:
        y = np.array([r["params"][nm] for r in ok], float)
        s = {"n": len(ok)}
        if nm in getattr(model, "proportional", ()) or nm in model.linear:
            # proporcional à potência: inclinação pela origem e reta completa
            a, ea = _line_origin(x, y)
            s.update(kind="proporcional", slope0=a, slope0_err=ea, line=_line(x, y),
                     ratio=(y / x).tolist())
        else:
            s.update(kind="constante", mean=float(y.mean()),
                     std=float(y.std(ddof=1)) if len(y) > 1 else float("nan"),
                     sem=float(y.std(ddof=1) / np.sqrt(len(y))) if len(y) > 1 else float("nan"),
                     line=_line(x, y))
            ln = s["line"]
            s["trend"] = bool(ln and abs(ln["slope"]) > 2 * ln["slope_err"])
        out[nm] = s
    return out


class BatchJob:
    """Série rodando em segundo plano no servidor (mesma interface de Job)."""

    def __init__(self, model: Model, params: dict, free_specs, items, options, acc, slots=None):
        self.id = uuid.uuid4().hex[:10]
        self.model, self.params, self.free_specs = model, params, free_specs
        self.items = [dict(name=str(it.get("name", f"transiente {k + 1}")),
                           t=np.asarray(it["t"], float), y=np.asarray(it["y"], float),
                           x=None if it.get("x") in (None, "") else float(it["x"]))
                      for k, it in enumerate(items)]
        self.options, self.acc, self.slots = dict(options or {}), acc, slots
        self.stop = threading.Event()
        self.lock = threading.Lock()
        self.status = "queued"
        self.started = time.time()
        self.log = []
        self.results = [dict(name=it["name"], x=it["x"], status="pendente") for it in self.items]
        self.stats = {}
        self.item, self.phase, self.iteration, self.n_iter = 0, "", 0, 0

    def say(self, msg):
        with self.lock:
            self.log.append(f"[{time.time() - self.started:6.1f} s] {msg}")

    def _progress(self, info):
        with self.lock:
            self.item = info.get("item", self.item)
            self.phase = info.get("phase", "")
            self.iteration = info.get("iteration", 0)
            self.n_iter = info.get("n_iter", 0)
            if self.results[self.item]["status"] == "pendente":
                self.results[self.item]["status"] = "ajustando"

    def _result(self, k, res):
        stats = series_stats(self.model, [r for r in self.results[:k]] + [res]
                             + self.results[k + 1:])
        with self.lock:
            self.results[k] = res
            self.stats = stats

    def run(self):
        if self.slots is not None:
            self.slots.acquire()
        try:
            with self.lock:
                self.status = "running"
                self.started = time.time()
            n = len(self.items)
            free = ", ".join(f["name"] for f in self.free_specs)
            self.say(f"Série com {n} transientes; parâmetros livres: {free}")
            fit_series(self.model, self.items, self.params, self.free_specs, self.options,
                       self.acc, callback=self._progress, stop=self.stop, log=self.say,
                       on_result=self._result)
            with self.lock:
                self.status = "done"
            self.say(f"Série concluída em {time.time() - self.started:.1f} s")
        except Stopped:
            self.say("Série interrompida; os transientes já ajustados foram mantidos.")
            with self.lock:
                self.status = "stopped"
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            self.say(f"ERRO: {e}")
            with self.lock:
                self.status = "error"
        finally:
            if self.slots is not None:
                self.slots.release()

    def snapshot(self):
        with self.lock:
            return dict(id=self.id, kind="serie", status=self.status, item=self.item,
                        n_items=len(self.items), phase=self.phase, iteration=self.iteration,
                        n_iter=self.n_iter, results=[dict(r) for r in self.results],
                        stats=self.stats, log=list(self.log),
                        elapsed=time.time() - self.started)
