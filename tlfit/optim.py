"""
Parameter fitting.

* Automatic mode (default): two-level fit. The θ (and the scale I0) enter the
  phase linearly, Φ = Σ θ_j Φ_j(D, …). The free nonlinear parameters (up to two)
  are swept on a grid; at each point Φ_j is computed once and the θ are fitted
  inside, from several starting points, at low cost. The best minima of the grid
  are refined by Levenberg–Marquardt with all parameters. Deterministic, and it
  covers the whole range of the nonlinear parameters.
* Particle swarm (PSO) + Levenberg–Marquardt, when there are more than two free
  nonlinear parameters (or when chosen in the interface).

The search uses a lighter numerical grid (Accuracy.search); the final refinement
and the uncertainties use the normal accuracy.

The free parameters are mapped to z ∈ [0, 1]. Positive parameters whose range
covers more than two decades (diffusivities, for example) use a logarithmic
scale, which lets the swarm explore orders of magnitude evenly.

Log messages are English ``Msg`` objects; the server translates them (tlfit/i18n.py).
"""
from __future__ import annotations

import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .i18n import Msg
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
        raise ValueError(Msg("equal bounds for {name}", name=name))
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
        """Counts evaluations and keeps the best point (for the display)."""
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
    """PSO with decreasing inertia (0.9 → 0.4) and global topology.

    callback(info) is called at every iteration with positions, χ² and best point."""
    rng = np.random.default_rng(seed)
    nd = len(prob.free)
    # initialization: Latin hypercube + the user's initial guess
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
            V[out] *= -0.5            # reflect at the walls
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
    """Local refinement (trust region with bounds) + uncertainties."""
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
    """Standard deviations (1σ) and correlation matrix of the free parameters."""
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
# Two-level fit
# ---------------------------------------------------------------------------

def split_free(prob: FitProblem):
    """(linear, nonlinear) among the free parameters."""
    lin_names = set(prob.model.linear) | {"scale"}
    lin = [i for i, fp in enumerate(prob.free) if fp.name in lin_names]
    nl = [i for i, fp in enumerate(prob.free) if fp.name not in lin_names]
    return lin, nl


def two_level_applicable(prob: FitProblem, max_nonlinear=2):
    if not getattr(prob.model, "phase_based", True):
        return False
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
    """Indices of the local minima of the grid (immediate neighbours), from best to worst."""
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
    """Optimal θ (and scale) for fixed nonlinear parameters. Returns (χ², full z, curve)."""
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
        # several starting points for the θ (the phase is periodic in them: there are local
        # minima); the scale I0 is well behaved and always starts from the current value.
        # The starting points are chosen by the phase they produce (±1, ±4 and ±15 rad at
        # most), not uniformly within the bounds: with wide bounds (θ_TM in 1/m), uniform
        # starts fell on phases of thousands of rad, physically meaningless and very
        # expensive to integrate.
        rng = np.random.default_rng(seed)
        n_starts = min(5, 1 + 2 * len(pos_phase))
        unit = pc.unit_phase()
        starts = [x0]
        for _ in range(n_starts - 1):
            xs = x0.copy()
            for j, k in enumerate(pos_phase):
                mag = rng.choice([1.0, 4.0, 15.0]) / max(unit[j], 1e-300)
                xs[k] = np.clip(rng.choice([-1.0, 1.0]) * mag, lo[k], hi[k])
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


SWEEP_POINTS = 200


def _sweep_problem(prob: FitProblem):
    """For long transients the sweep uses ~200 points equally spaced in index
    (the final fit uses all of them): the cost of the grid drops in proportion. Evaluations
    are counted in the original problem; the curve is only kept in stages with all points."""
    n = len(prob.t)
    if n <= SWEEP_POINTS * 1.25:
        return prob
    idx = np.unique(np.round(np.linspace(0, n - 1, SWEEP_POINTS)).astype(int))
    sub = FitProblem(prob.model, prob.t[idx], prob.y[idx], prob.params, prob.free, prob.acc)
    factor = n / len(idx)

    def record(z, chi, f, n=1):
        prob.record(z, chi * factor, None, n)
    sub.record = record
    return sub


def two_level(prob: FitProblem, callback=None, stop=None, log=print, workers=None):
    """Sweep over the nonlinear parameters + θ inside; returns sorted candidates (z, χ²)."""
    lin, nl = split_free(prob)
    Znl = _grid_points(len(nl))
    acc = prob.search_acc
    names_nl = ", ".join(prob.free[i].name for i in nl) or "—"
    names_lin = ", ".join(prob.free[i].name for i in lin) or "—"
    log(Msg("Two-level fit: sweep over {nl} ({n} points); {lin} fitted at each point",
            nl=names_nl, n=len(Znl), lin=names_lin))
    N = len(Znl)
    F = np.full(N, np.inf)
    Zfull = np.zeros((N, len(prob.free)))
    done = []
    workers = workers or min(N, os.cpu_count() or 1)
    pool = ThreadPoolExecutor(max_workers=workers)
    lock = threading.Lock()
    sweep = _sweep_problem(prob)
    try:
        futs = {pool.submit(inner_fit, sweep, Znl[k], lin, nl, acc, 1000 + k, stop): k
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
                    callback(dict(phase="grid", iteration=len(done), n_iter=N,
                                  positions=Zfull[idx].copy(), chi=F[idx].copy(),
                                  gbest=prob.best_z, gchi=prob.best_chi))
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
    if sweep is not prob:
        # the sweep χ² is estimated on few points: it does not count as the "best" of the full fit
        prob.best_chi, prob.best_z, prob.best_curve = np.inf, None, None
    mins = _local_minima(Znl, F, len(nl))
    best = F[mins[0]] if mins else np.inf
    cand = [i for i in mins if F[i] <= 3 * best][:3]
    if len(mins) > 1:
        log(Msg("Sweep: {n} local minima; refining {k} (χ² {chis})", n=len(mins), k=len(cand),
                chis=", ".join(f"{F[i]:.3e}" for i in cand)))
    else:
        log(Msg("Sweep: a single minimum (χ² {chi:.3e})", chi=best))
    return [(Zfull[i], F[i]) for i in cand]


def run_fit(prob: FitProblem, method="auto", n_particles=20, n_iter=40, seed=None,
            callback=None, stop=None, log=print, restarts=3):
    """Runs the complete fit and returns a dictionary with the result.

    method: "auto" (two levels when possible, otherwise PSO+LM), "pso+lm", "pso", "lm"."""
    t0 = time.time()
    z0 = np.array([fp.to_z(prob.params[fp.name]) for fp in prob.free])
    if method == "auto":
        if two_level_applicable(prob):
            method = "grid"
        else:
            _, nl = split_free(prob)
            log(Msg("{n} free nonlinear parameters ({names}): using swarm (PSO) + LM",
                    n=len(nl), names=", ".join(prob.free[i].name for i in nl)))
            method = "pso+lm"
    use_pso = method in ("pso", "pso+lm")
    use_lm = method in ("lm", "pso+lm", "grid")
    rng = np.random.default_rng(seed)
    runs = []   # (χ², z, sigma, extra)

    if method == "grid":
        cands = two_level(prob, callback=callback, stop=stop, log=log)

        def cb(info):
            if callback:
                callback(dict(info, run=0, n_runs=1))
        for zc, _ in cands:
            res = levenberg_marquardt(prob, zc, callback=cb, stop=stop)
            sig, extra = uncertainties(prob, res)
            runs.append((float(res.fun @ res.fun), res.x, sig, extra))
        log(Msg("Levenberg–Marquardt (normal accuracy, all parameters): χ² = {chi:.4e}",
                chi=min(r[0] for r in runs)))
    else:
        n_runs = max(1, int(restarts)) if use_pso else 1
        for run in range(n_runs):
            tag = Msg("Run {run} of {n}: ", run=run + 1, n=n_runs) if n_runs > 1 else ""
            z = z0

            def cb(info, run=run):
                if callback:
                    callback(dict(info, run=run, n_runs=n_runs))

            if use_pso:
                log(Msg("{tag}PSO with {particles} particles × {iters} iterations",
                        tag=tag, particles=n_particles, iters=n_iter))
                z, chi = pso(prob, n_particles, n_iter, seed=int(rng.integers(2 ** 31)),
                             callback=cb, stop=stop, x0_z=z0)
                log(Msg("{tag}PSO finished with χ² = {chi:.4e}", tag=tag, chi=chi))
            sig, extra = {}, None
            if use_lm:
                res = levenberg_marquardt(prob, z, callback=cb, stop=stop)
                z = res.x
                sig, extra = uncertainties(prob, res)
                log(Msg("{tag}Levenberg–Marquardt finished with χ² = {chi:.4e}",
                        tag=tag, chi=float(res.fun @ res.fun)))
            runs.append((prob.chi2(z), z, sig, extra))

        chis = np.array([r[0] for r in runs])
        if n_runs > 1:
            best = int(np.argmin(chis))
            same = int(np.sum(chis <= chis[best] * (1 + 1e-3) + 1e-300))
            if same == n_runs:
                log(Msg("All {n} runs reached the same minimum.", n=n_runs))
            else:
                log(Msg("Warning: only {same} of {n} runs reached the best minimum (the others "
                        "stopped at χ² = {chis}). There is more than one minimum; if the result "
                        "looks odd, increase the number of runs or narrow the bounds.",
                        same=same, n=n_runs,
                        chis=", ".join(f"{c:.3e}" for c in np.delete(chis, best))))

    chi_best, z, sig, extra = min(runs, key=lambda r: r[0])
    params = prob.phys(z)
    curve = prob.curve(z)
    chi = float(((curve - prob.y) ** 2).sum())
    with prob._lock:   # the display shows the final result (normal accuracy)
        prob.best_chi, prob.best_z, prob.best_curve = chi, np.array(z, float), curve
    log(Msg("Total: {n} model evaluations in {s:.1f} s", n=prob.n_eval, s=time.time() - t0))
    return dict(params=params, chi2=chi, sigma=sig, extra=extra, curve=curve, z=z,
                n_eval=prob.n_eval, seconds=time.time() - t0)
