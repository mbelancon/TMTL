"""
Ajuste de parâmetros.

* Modo automático (padrão): ajuste em dois níveis. Os θ (e a escala I0) entram
  linearmente na fase, Φ = Σ θ_j Φ_j(D, …). Varre-se os parâmetros não lineares
  livres (até dois) numa grade; em cada ponto, Φ_j é calculada uma vez e os θ são
  ajustados por dentro, com vários pontos de partida, a um custo baixo. Os melhores
  mínimos da grade são refinados por Levenberg–Marquardt com todos os parâmetros.
  Determinístico e cobre toda a faixa dos parâmetros não lineares.
* Enxame de partículas (PSO) + Levenberg–Marquardt, quando há mais de dois
  parâmetros não lineares livres (ou se escolhido na interface).

A busca usa uma malha numérica mais leve (Accuracy.search); o refinamento final
e as incertezas usam a precisão normal.

Os parâmetros livres são mapeados em z ∈ [0, 1]. Parâmetros positivos cujo
intervalo cobre mais de duas décadas (difusividades, por exemplo) usam escala
logarítmica, o que deixa o enxame explorar ordens de grandeza por igual.
"""
from __future__ import annotations

import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .models import Accuracy, Model, PhaseCache


@dataclass
class FreeParam:
    name: str
    lo: float
    hi: float
    log: bool = False

    def to_phys(self, z):
        z = np.clip(z, 0.0, 1.0)
        if self.log:
            return self.lo * (self.hi / self.lo) ** z
        return self.lo + (self.hi - self.lo) * z

    def to_z(self, x):
        if self.log:
            return float(np.clip(np.log(x / self.lo) / np.log(self.hi / self.lo), 0, 1))
        return float(np.clip((x - self.lo) / (self.hi - self.lo), 0, 1))

    def dphys_dz(self, x):
        return x * np.log(self.hi / self.lo) if self.log else (self.hi - self.lo)


def make_free(name, lo, hi, log=None):
    lo, hi = float(lo), float(hi)
    if hi < lo:
        lo, hi = hi, lo
    if hi == lo:
        raise ValueError(f"limites iguais para {name}")
    if log is None:
        log = lo > 0 and hi / lo > 100
    return FreeParam(name, lo, hi, bool(log and lo > 0))


@dataclass
class FitProblem:
    model: Model
    t: np.ndarray
    y: np.ndarray
    params: dict
    free: list[FreeParam]
    acc: Accuracy = field(default_factory=Accuracy)

    def __post_init__(self):
        self.t = np.ascontiguousarray(self.t, float)
        self.y = np.ascontiguousarray(self.y, float)
        self.search_acc = self.acc.search()
        self._lock = threading.Lock()
        self.n_eval = 0
        self.best_chi = np.inf
        self.best_z = None
        self.best_curve = None

    def phys(self, z) -> dict:
        p = dict(self.params)
        for fp, zi in zip(self.free, z):
            p[fp.name] = float(fp.to_phys(zi))
        return p

    def curve(self, z, acc=None):
        return self.model.intensity(self.t, self.phys(z), acc or self.acc)

    def record(self, z, chi, f, n=1):
        """Contabiliza avaliações e guarda o melhor ponto (para a tela)."""
        with self._lock:
            self.n_eval += n
            if chi < self.best_chi:
                self.best_chi, self.best_z, self.best_curve = chi, np.array(z, float), f

    def residuals(self, z, acc=None):
        f = self.curve(z, acc)
        r = f - self.y
        if not np.all(np.isfinite(r)):
            r = np.where(np.isfinite(r), r, 1e3)
        self.record(z, float(r @ r), f)
        return r

    def chi2(self, z):
        r = self.residuals(z)
        return float(r @ r)

    def chi2_search(self, z):
        r = self.residuals(z, self.search_acc)
        return float(r @ r)


class Stopped(Exception):
    pass


def pso(prob: FitProblem, n_particles=20, n_iter=40, seed=None, workers=None,
        callback=None, stop: threading.Event | None = None, x0_z=None):
    """PSO com inércia decrescente (0.9 → 0.4) e topologia global.

    callback(info) é chamado a cada iteração com posições, χ² e melhor ponto."""
    rng = np.random.default_rng(seed)
    nd = len(prob.free)
    # inicialização: hipercubo latino + a estimativa inicial do usuário
    X = (rng.permuted(np.tile(np.arange(n_particles), (nd, 1)), axis=1).T
         + rng.random((n_particles, nd))) / n_particles
    if x0_z is not None:
        X[0] = x0_z
    V = (rng.random((n_particles, nd)) - 0.5) * 0.2
    vmax = 0.25
    c1 = c2 = 1.49445
    workers = workers or min(n_particles, os.cpu_count() or 1)
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        def evaluate(P):
            return np.array(list(pool.map(prob.chi2_search, P)))

        F = evaluate(X)
        pbest, pf = X.copy(), F.copy()
        g = int(np.argmin(pf))
        gbest, gf = pbest[g].copy(), pf[g]
        for it in range(n_iter):
            if stop is not None and stop.is_set():
                raise Stopped
            if callback:
                callback(dict(phase="pso", iteration=it, n_iter=n_iter, positions=X.copy(),
                              chi=F.copy(), gbest=gbest.copy(), gchi=gf))
            w = 0.9 - 0.5 * it / max(1, n_iter - 1)
            r1, r2 = rng.random((2, n_particles, nd))
            V = w * V + c1 * r1 * (pbest - X) + c2 * r2 * (gbest - X)
            V = np.clip(V, -vmax, vmax)
            X = X + V
            out = (X < 0) | (X > 1)
            X = np.clip(X, 0, 1)
            V[out] *= -0.5            # reflete nas paredes
            F = evaluate(X)
            better = F < pf
            pbest[better], pf[better] = X[better], F[better]
            g = int(np.argmin(pf))
            if pf[g] < gf:
                gbest, gf = pbest[g].copy(), pf[g]
        if callback:
            callback(dict(phase="pso", iteration=n_iter, n_iter=n_iter, positions=X.copy(),
                          chi=F.copy(), gbest=gbest.copy(), gchi=gf))
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    return gbest, gf


def levenberg_marquardt(prob: FitProblem, z0, max_nfev=300, callback=None,
                        stop: threading.Event | None = None):
    """Refinamento local (região de confiança com limites) + incertezas."""
    def fun(z):
        if stop is not None and stop.is_set():
            raise Stopped
        r = prob.residuals(z)
        if callback:
            callback(dict(phase="lm", n_eval=prob.n_eval, gchi=prob.best_chi,
                          gbest=prob.best_z))
        return r

    z0 = np.clip(np.asarray(z0, float), 1e-9, 1 - 1e-9)
    res = least_squares(fun, z0, bounds=(0.0, 1.0), method="trf", x_scale="jac",
                        diff_step=1e-6, ftol=1e-12, xtol=1e-12, gtol=1e-12,
                        max_nfev=max_nfev)
    return res


def uncertainties(prob: FitProblem, res):
    """Desvios-padrão (1σ) e matriz de correlação dos parâmetros livres."""
    z = res.x
    J = res.jac
    n, p = J.shape
    chi = float(res.fun @ res.fun)
    dof = max(1, n - p)
    s2 = chi / dof
    names = [fp.name for fp in prob.free]
    try:
        cov_z = np.linalg.pinv(J.T @ J) * s2
    except np.linalg.LinAlgError:
        return {nm: float("nan") for nm in names}, None
    phys = prob.phys(z)
    D = np.array([fp.dphys_dz(phys[fp.name]) for fp in prob.free])
    cov = cov_z * np.outer(D, D)
    sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(sig, sig)
    at_bound = [nm for nm, zi in zip(names, z) if zi < 1e-6 or zi > 1 - 1e-6]
    return dict(zip(names, sig.tolist())), dict(corr=np.nan_to_num(corr).tolist(),
                                                 at_bound=at_bound, reduced_chi2=s2, dof=dof)


# ---------------------------------------------------------------------------
# Ajuste em dois níveis
# ---------------------------------------------------------------------------

def split_free(prob: FitProblem):
    """(lineares, não lineares) entre os parâmetros livres."""
    lin_names = set(prob.model.linear) | {"scale"}
    lin = [i for i, fp in enumerate(prob.free) if fp.name in lin_names]
    nl = [i for i, fp in enumerate(prob.free) if fp.name not in lin_names]
    return lin, nl


def two_level_applicable(prob: FitProblem, max_nonlinear=2):
    _, nl = split_free(prob)
    return len(nl) <= max_nonlinear


def _grid_points(n_nl):
    if n_nl == 0:
        return np.zeros((1, 0))
    if n_nl == 1:
        return np.linspace(0, 1, 40)[:, None]
    g = np.linspace(0, 1, 12)
    return np.array([[a, b] for a in g for b in g])


def _local_minima(Z, F, n_nl):
    """Índices dos mínimos locais da grade (vizinhos imediatos), do melhor ao pior."""
    order = np.argsort(F)
    if n_nl == 0:
        return [0]
    step = (Z[1, -1] - Z[0, -1]) if n_nl == 1 else (1 / 11)
    mins = []
    for i in order:
        if not np.isfinite(F[i]):
            break
        near = np.all(np.abs(Z - Z[i]) <= step * 1.01, axis=1)
        if F[i] <= F[near].min():
            mins.append(int(i))
    return mins


def inner_fit(prob: FitProblem, z_nl, lin, nl, acc, seed, stop=None):
    """θ ótimos (e escala) para parâmetros não lineares fixos. Devolve (χ², z completo, curva)."""
    free = prob.free
    p = dict(prob.params)
    for i, zi in zip(nl, z_nl):
        p[free[i].name] = float(free[i].to_phys(zi))
    phase_lin = [i for i in lin if free[i].name != "scale"]
    scale_i = next((i for i in lin if free[i].name == "scale"), None)
    pc = PhaseCache(prob.model, prob.t, p, [free[i].name for i in phase_lin], acc)
    z = np.zeros(len(free))
    z[nl] = z_nl
    if not pc.ok:
        prob.record(z, np.inf, None)
        return np.inf, z, None
    lo = np.array([free[i].lo for i in lin])
    hi = np.array([free[i].hi for i in lin])
    x0 = np.clip([p[free[i].name] for i in lin], lo, hi)
    n_phase = len(phase_lin)
    fixed_scale = p.get("scale", 1.0)

    pos_scale = lin.index(scale_i) if scale_i is not None else None
    pos_phase = [lin.index(i) for i in phase_lin]

    def unpack(x):
        th = [x[k] for k in pos_phase]
        sc = fixed_scale if pos_scale is None else x[pos_scale]
        return th, sc

    def model_curve(x):
        th, sc = unpack(x)
        return pc.intensity(th, sc)

    nfev = 0
    cache = {}

    def evaluate(x):
        nonlocal nfev
        key = x.tobytes()
        if cache.get("key") != key:
            if stop is not None and stop.is_set():
                raise Stopped
            nfev += 1
            th, sc = unpack(x)
            I, J = pc.intensity(th, sc, jac=True)
            r = I - prob.y
            bad = ~np.isfinite(r)
            r[bad] = 1e3
            if J is not None:
                Jx = np.zeros((len(r), len(lin)))
                for k, pk in enumerate(pos_phase):
                    Jx[:, pk] = J[:, k]
                if pos_scale is not None:
                    Jx[:, pos_scale] = J[:, -1]
                Jx[bad] = 0
            else:
                Jx = np.zeros((len(r), len(lin)))
            cache.update(key=key, r=r, J=Jx)
        return cache

    def res(x):
        return evaluate(np.asarray(x, float))["r"]

    def jac(x):
        return evaluate(np.asarray(x, float))["J"]

    if not lin:
        x_best = x0
        r = res(x0)
        chi = float(r @ r)
    else:
        # vários pontos de partida para os θ (a fase é periódica neles: há mínimos
        # locais); a escala I0 é bem comportada e parte sempre do valor atual
        rng = np.random.default_rng(seed)
        n_starts = min(5, 1 + 2 * len(pos_phase))
        starts = [x0]
        for _ in range(n_starts - 1):
            xs = x0.copy()
            for k in pos_phase:
                xs[k] = lo[k] + (hi[k] - lo[k]) * rng.random()
            starts.append(xs)
        chi, x_best = np.inf, x0
        for xs in starts:
            r = least_squares(res, xs, jac=jac, bounds=(lo, hi), x_scale=(hi - lo) / 10,
                              method="trf", ftol=1e-10, xtol=1e-10, max_nfev=80)
            c = float(r.fun @ r.fun)
            if c < chi:
                chi, x_best = c, r.x
    for k, i in enumerate(lin):
        z[i] = free[i].to_z(x_best[k])
    f = model_curve(x_best)
    prob.record(z, chi, f, n=nfev)
    return chi, z, f


def two_level(prob: FitProblem, callback=None, stop=None, log=print, workers=None):
    """Varredura nos não lineares + θ por dentro; devolve candidatos (z, χ²) ordenados."""
    lin, nl = split_free(prob)
    Znl = _grid_points(len(nl))
    acc = prob.search_acc
    names_nl = ", ".join(prob.free[i].name for i in nl) or "nenhum"
    names_lin = ", ".join(prob.free[i].name for i in lin) or "nenhum"
    log(f"Ajuste em dois níveis: varredura em {names_nl} ({len(Znl)} pontos); "
        f"{names_lin} ajustados em cada ponto")
    N = len(Znl)
    F = np.full(N, np.inf)
    Zfull = np.zeros((N, len(prob.free)))
    done = []
    workers = workers or min(N, os.cpu_count() or 1)
    pool = ThreadPoolExecutor(max_workers=workers)
    lock = threading.Lock()
    try:
        futs = {pool.submit(inner_fit, prob, Znl[k], lin, nl, acc, 1000 + k, stop): k
                for k in range(N)}
        from concurrent.futures import as_completed
        for fut in as_completed(futs):
            k = futs[fut]
            chi, z, _ = fut.result()
            with lock:
                F[k], Zfull[k] = chi, z
                done.append(k)
                if callback:
                    idx = np.array(done)
                    callback(dict(phase="grade", iteration=len(done), n_iter=N,
                                  positions=Zfull[idx].copy(), chi=F[idx].copy(),
                                  gbest=prob.best_z, gchi=prob.best_chi))
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    mins = _local_minima(Znl, F, len(nl))
    best = F[mins[0]] if mins else np.inf
    cand = [i for i in mins if F[i] <= 3 * best][:3]
    if len(mins) > 1:
        log(f"Varredura: {len(mins)} mínimos locais; refinando {len(cand)} "
            f"(χ² {', '.join(f'{F[i]:.3e}' for i in cand)})")
    else:
        log(f"Varredura: um único mínimo (χ² {best:.3e})")
    return [(Zfull[i], F[i]) for i in cand]


def run_fit(prob: FitProblem, method="auto", n_particles=20, n_iter=40, seed=None,
            callback=None, stop=None, log=print, restarts=3):
    """Executa o ajuste completo e devolve um dicionário com o resultado.

    method: "auto" (dois níveis quando possível, senão PSO+LM), "pso+lm", "pso", "lm"."""
    t0 = time.time()
    z0 = np.array([fp.to_z(prob.params[fp.name]) for fp in prob.free])
    if method == "auto":
        if two_level_applicable(prob):
            method = "grade"
        else:
            _, nl = split_free(prob)
            log(f"{len(nl)} parâmetros não lineares livres "
                f"({', '.join(prob.free[i].name for i in nl)}): usando enxame (PSO) + LM")
            method = "pso+lm"
    use_pso = method in ("pso", "pso+lm")
    use_lm = method in ("lm", "pso+lm", "grade")
    rng = np.random.default_rng(seed)
    runs = []   # (χ², z, sigma, extra)

    if method == "grade":
        cands = two_level(prob, callback=callback, stop=stop, log=log)

        def cb(info):
            if callback:
                callback(dict(info, run=0, n_runs=1))
        for zc, _ in cands:
            res = levenberg_marquardt(prob, zc, callback=cb, stop=stop)
            sig, extra = uncertainties(prob, res)
            runs.append((float(res.fun @ res.fun), res.x, sig, extra))
        log(f"Levenberg–Marquardt (precisão normal, todos os parâmetros): χ² = "
            f"{min(r[0] for r in runs):.4e}")
    else:
        n_runs = max(1, int(restarts)) if use_pso else 1
        for run in range(n_runs):
            tag = f"Rodada {run + 1} de {n_runs}: " if n_runs > 1 else ""
            z = z0

            def cb(info, run=run):
                if callback:
                    callback(dict(info, run=run, n_runs=n_runs))

            if use_pso:
                log(f"{tag}PSO com {n_particles} partículas × {n_iter} iterações")
                z, chi = pso(prob, n_particles, n_iter, seed=int(rng.integers(2 ** 31)),
                             callback=cb, stop=stop, x0_z=z0)
                log(f"{tag}PSO terminou com χ² = {chi:.4e}")
            sig, extra = {}, None
            if use_lm:
                res = levenberg_marquardt(prob, z, callback=cb, stop=stop)
                z = res.x
                sig, extra = uncertainties(prob, res)
                log(f"{tag}Levenberg–Marquardt terminou com χ² = {float(res.fun @ res.fun):.4e}")
            runs.append((prob.chi2(z), z, sig, extra))

        chis = np.array([r[0] for r in runs])
        if n_runs > 1:
            best = int(np.argmin(chis))
            same = int(np.sum(chis <= chis[best] * (1 + 1e-3) + 1e-300))
            if same == n_runs:
                log(f"As {n_runs} rodadas chegaram ao mesmo mínimo.")
            else:
                log(f"Atenção: só {same} de {n_runs} rodadas chegaram ao melhor mínimo "
                    f"(as outras pararam em χ² = "
                    f"{', '.join(f'{c:.3e}' for c in np.delete(chis, best))}). "
                    "Há mais de um mínimo; se o resultado parecer estranho, aumente o número "
                    "de rodadas ou estreite os limites.")

    chi_best, z, sig, extra = min(runs, key=lambda r: r[0])
    params = prob.phys(z)
    curve = prob.curve(z)
    chi = float(((curve - prob.y) ** 2).sum())
    with prob._lock:   # a tela mostra o resultado final (precisão normal)
        prob.best_chi, prob.best_z, prob.best_curve = chi, np.array(z, float), curve
    log(f"Total: {prob.n_eval} avaliações do modelo em {time.time() - t0:.1f} s")
    return dict(params=params, chi2=chi, sigma=sig, extra=extra, curve=curve, z=z,
                n_eval=prob.n_eval, seconds=time.time() - t0)
