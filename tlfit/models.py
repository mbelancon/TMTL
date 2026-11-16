"""
Time-resolved thermal lens / thermal mirror models.

Vectorized reimplementation of the TMTL models (tllam.c and tmhamtl-onoff.c).
The integrand formulas are exactly those of the original C code; what changes
is how the integrals are computed:

  I(t) = (1+v^2) | ∫_0^Gmax e^{-g} exp(i[v g + Φ(g,t)]) dg |^2
  Φ(g,t) = ∫ K(α,t) [J0(α w0 sqrt(m g)) - 1] dα

The integrand of Φ is separable: K(α,t) depends on time and on the thermal
parameters; the Bessel factor depends only on (α, g, m, w0). With fixed
quadratures (panel Gauss-Legendre) in α and in g, Φ becomes a matrix product

  Φ[t, g] = K[t, α] @ M[α, g],   M = w_α (J0 - 1)

and M is computed only once per geometry (m, w0). The real and imaginary
parts of the amplitude come out of the same complex sum. In the original code
each Φ was an adaptive integral recomputed inside another adaptive integral
(and twice: once for cos, once for sin).

The mathematics of every model is described in docs/MATH_METHODS.md.
User-facing strings (titles, descriptions, labels, θ definitions) are written
in English; their Portuguese translations live in tlfit/i18n.py.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy.special import erf, erfc, erfcx, dawsn, j0, expit

SQRT_PI = np.sqrt(np.pi)

# ---------------------------------------------------------------------------
# Quadratures
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _gl(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def composite_gl(edges, n):
    """n-point Gauss-Legendre on each panel [edges[i], edges[i+1]]."""
    x, w = _gl(n)
    edges = np.asarray(edges, float)
    a, b = edges[:-1, None], edges[1:, None]
    nodes = 0.5 * (b - a) * x[None, :] + 0.5 * (b + a)
    weights = 0.5 * (b - a) * w[None, :]
    return nodes.ravel(), weights.ravel()


@dataclass(frozen=True)
class Accuracy:
    """Resolution control of the quadratures (factor 1 = default)."""
    u_max: float = 14.0        # cut at u = α w0 (e^{-u²/8} ≈ 2e-11)
    g_max: float = 30.0        # same limit as the original code
    nodes_per_panel: int = 12
    u_panels_per_period: float = 1.0
    g_panels_per_period: float = 1.0
    u_grading_levels: int = 14  # geometric panels near u = 0
    tau_nodes: int = 80         # tanh-sinh for the τ integral (air, dn/dT)

    def scaled(self, f: float) -> "Accuracy":
        return Accuracy(self.u_max, self.g_max, self.nodes_per_panel,
                        self.u_panels_per_period * f, self.g_panels_per_period * f,
                        self.u_grading_levels, int(self.tau_nodes * f))

    def search(self) -> "Accuracy":
        """Lighter version, for the search phase of the fit.

        g up to 18 (the neglected term is ~2e^{-18} ≈ 3e-8 relative) and 0.7× grids:
        error in I(t) ≲ 3e-7, far below typical noise (~1e-4), with 2–3× less
        work. The final refinement always uses the normal accuracy."""
        return Accuracy(self.u_max, min(self.g_max, 18.0), self.nodes_per_panel,
                        self.u_panels_per_period * 0.7, self.g_panels_per_period * 0.7,
                        self.u_grading_levels, max(40, int(self.tau_nodes * 0.7)))


def u_grid(m: float, w0: float, acc: Accuracy):
    u_top = min(acc.u_max, 2e5 * w0)  # the original integrates α over [0, 2e5]
    # period of J0(u sqrt(m g)) in u in the worst case (g = g_max)
    period = 2 * np.pi / np.sqrt(max(m, 1e-12) * acc.g_max)
    width = period / acc.u_panels_per_period
    first = min(width, u_top)
    graded = first * 2.0 ** -np.arange(acc.u_grading_levels, 0, -1)
    uniform = np.linspace(first, u_top, max(1, int(np.ceil((u_top - first) / width))) + 1)
    edges = np.concatenate([[0.0], graded, uniform])
    edges = np.unique(edges)
    return composite_gl(edges, acc.nodes_per_panel)


def g_base_edges(v: float, acc: Accuracy):
    """Base panels in g: width = one period of e^{i v g} (or 1, if v < 2π)."""
    rate = max(abs(v), 2 * np.pi)
    width = 2 * np.pi / rate / acc.g_panels_per_period
    n = int(np.ceil(acc.g_max / width))
    return np.linspace(0, acc.g_max, n + 1)


def g_grid(v: float, acc: Accuracy, split: tuple):
    """Base panels subdivided (split[i] = 1, 2, 4, ...) where the phase turns fast."""
    base = g_base_edges(v, acc)
    edges = [base[0]]
    for i, si in enumerate(split):
        edges.extend(np.linspace(base[i], base[i + 1], si + 1)[1:])
    return composite_gl(np.array(edges), acc.nodes_per_panel)


@lru_cache(maxsize=4)
def _alpha_grid(m, w0, acc: Accuracy):
    u, wu = u_grid(m, w0, acc)
    return u, wu


@lru_cache(maxsize=8)
def _bessel_matrix(m, w0, acc: Accuracy, v, split: tuple):
    u, wu = _alpha_grid(m, w0, acc)
    g, wg = g_grid(v, acc, split)
    M = (wu / w0)[:, None] * (j0(u[:, None] * np.sqrt(m * g[None, :])) - 1.0)
    return M, g, wg * np.exp(-g)


MAX_SPLIT = 64  # cap of sub-panels per base panel (phase > ~400 rad: unphysical)
_SPLIT_CACHE: dict = {}   # (m, w0, v, acc) → subdivision used in the last evaluation


def needed_split(phi, g, base_edges, n_panels):
    """Sub-panels needed in each base panel so that the total phase
    (v g + Φ) turns at most ~2π per panel. Powers of 2 (reuses the cache)."""
    with np.errstate(invalid="ignore"):
        rate = np.abs(np.diff(phi, axis=1)).max(axis=0) / np.diff(g)   # between neighbouring nodes
    mid = 0.5 * (g[1:] + g[:-1])
    idx = np.clip(np.searchsorted(base_edges, mid) - 1, 0, n_panels - 1)
    rmax = np.zeros(n_panels)
    np.maximum.at(rmax, idx, np.nan_to_num(rate, nan=np.inf))
    need = rmax * np.diff(base_edges) / (2 * np.pi)
    need = np.nan_to_num(need, nan=MAX_SPLIT, posinf=MAX_SPLIT)
    new = np.ones(n_panels, int)
    big = need > 1
    new[big] = 2 ** np.ceil(np.log2(np.minimum(need[big], MAX_SPLIT))).astype(int)
    return tuple(int(x) for x in new)


def solve_phase(phase_fn, m, w0, v, acc: Accuracy):
    """Chooses the g grid and returns (φ, g, weights).

    phase_fn(M, g) returns the total phase (Nt, Ng) on the given grid. It starts
    from the subdivision used in the previous evaluation with the same geometry
    (it almost always works, and then a single K·M product is enough) and
    refines where needed."""
    key = (m, w0, v, acc)
    base = g_base_edges(v, acc)
    n = len(base) - 1
    split = _SPLIT_CACHE.get(key, (1,) * n)
    for _ in range(8):
        M, g, wg = _bessel_matrix(m, w0, acc, v, split)
        phi = phase_fn(M, g)
        need = needed_split(phi, g, base, n)
        if all(a <= b for a, b in zip(need, split)):
            _SPLIT_CACHE[key] = need      # may shrink in the next evaluation
            return phi, g, wg
        split = tuple(max(a, b) for a, b in zip(need, split))
    return phi, g, wg


def refine_split(phi, g, v, base_edges, split, nodes_per_panel):
    """Compatibility: subdivision needed, never smaller than the current one."""
    need = needed_split(phi, g, base_edges, len(split))
    return tuple(max(a, b) for a, b in zip(need, split))


def tanh_sinh(n: int):
    """tanh-sinh nodes/weights on [0, 1] (good for boundary layers at the ends)."""
    h = 6.0 / n
    k = np.arange(-n // 2, n // 2 + 1) * h
    s = 0.5 * np.pi * np.sinh(k)
    x = expit(2 * s)
    w = 0.5 * h * 0.5 * np.pi * np.cosh(k) / np.cosh(s) ** 2
    keep = (x > 0) & (x < 1) & (w > 1e-300)
    return x[keep], w[keep]


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


@dataclass
class Param:
    name: str
    value: float
    label: str
    unit: str = ""
    fit: bool = False
    lo: float | None = None
    hi: float | None = None


@dataclass
class Model:
    key: str
    title: str
    description: str
    params: list[Param] = field(default_factory=list)
    # presets: name → {"values": {param: value}, "bounds": {param: [min, max]}}
    presets: dict = field(default_factory=dict)

    def defaults(self) -> dict[str, float]:
        return {p.name: p.value for p in self.params}

    # parameters that enter the phase linearly: Φ = Σ θ_j Φ_j (see PhaseCache).
    # Class attribute (no annotation, so it does not become a dataclass field).
    linear = ()
    # False for closed-form models (no K·M): the two-level fit does not apply
    phase_based = True
    # physical definition of each fittable θ: {param: {"def", "unit", "note"}} (shown in the interface)
    theta_info = {}

    # linear parameters whose phase does not have the K·M form (e.g. population lens): the
    # subclass returns Φ_j (Nt, Ng) for θ_j = 1 in direct_phase(name, t, g, p)
    direct = ()

    # implemented by the subclasses
    def kernel(self, t, alpha, p, acc):
        raise NotImplementedError

    def direct_phase(self, name, t, g, p):
        raise NotImplementedError

    def intensity(self, t, p: dict, acc: Accuracy = Accuracy(), info: dict | None = None):
        t = np.asarray(t, float)
        out = np.ones_like(t)
        live = t > 1e-9
        scale = p.get("scale", 1.0)
        if not np.any(live):
            return out * scale
        m, w0, v = float(p["m"]), float(p["w0"]), float(p["v"])
        u, _ = _alpha_grid(m, w0, acc)
        tl = t[live]
        with np.errstate(all="ignore"):
            K = self.kernel(tl[:, None], (u / w0)[None, :], p, acc)    # (Nt, Nα)
        if not np.all(np.isfinite(K)):
            # parameters outside the model's domain: return NaN (the fit discards it)
            out[live] = np.nan
            if info is not None:
                info.update(n_alpha=len(u), n_g=0, max_phase=float("nan"))
            return out * scale
        def phase_fn(M, g):
            ph = K @ M + v * g[None, :]
            for n in self.direct:
                if p.get(n, 0.0) != 0:
                    ph = ph + p[n] * self.direct_phase(n, tl[:, None], g, p)
            return ph
        phi, g, wg = solve_phase(phase_fn, m, w0, v, acc)
        amp = np.cos(phi) @ wg + 1j * (np.sin(phi) @ wg)
        out[live] = (amp.real ** 2 + amp.imag ** 2) * (1 + v * v)
        if info is not None:
            info.update(n_alpha=len(u), n_g=len(g), max_phase=float(np.abs(phi - v * g).max()))
        return out * scale


class PhaseCache:
    """Phase decomposed in the linear parameters, for fixed nonlinear parameters.

    Φ(t, g) = Φ_0 + Σ_j θ_j Φ_j, with Φ_j = K_j · M and K_j the kernel with θ_j = 1 and
    the other linear θ = 0. With this, varying the θ (and the scale I0) costs only
    the cos/sin and the sum in g, 10–15× less than a full model evaluation."""

    def __init__(self, model: Model, t, p: dict, lin_names, acc: Accuracy):
        self.t = np.asarray(t, float)
        self.live = self.t > 1e-9
        self.lin = list(lin_names)
        self.p = dict(p)
        self.acc = acc
        self.model = model
        self.m, self.w0, self.v = float(p["m"]), float(p["w0"]), float(p["v"])
        u, _ = _alpha_grid(self.m, self.w0, acc)
        tl = self.t[self.live][:, None]
        self.tl = tl
        a = (u / self.w0)[None, :]
        zero = {n: 0.0 for n in model.linear}
        with np.errstate(all="ignore"):
            # fixed part: linear θ that are NOT being varied, with their values
            p0 = {**p, **{n: 0.0 for n in self.lin}}
            self.K0 = model.kernel(tl, a, p0, acc) if any(p0[n] != 0 for n in model.linear) \
                else np.zeros((tl.shape[0], a.shape[1]))
            self.Kj = [model.kernel(tl, a, {**p, **zero, n: 1.0}, acc) for n in self.lin]
        self.ok = bool(np.all(np.isfinite(self.K0)) and all(np.all(np.isfinite(k)) for k in self.Kj))
        self._grid = None   # (grid id, Φ_0, [Φ_j], g, wg)

    def _set_grid(self, split):
        """Φ_0, Φ_j and the maximum rates of change of each Φ_j per base panel."""
        M, g, wg = _bessel_matrix(self.m, self.w0, self.acc, self.v, split)
        P0 = self.K0 @ M
        direct = self.model.direct
        for n in direct:                      # fixed direct terms (not varied)
            if n not in self.lin and self.p.get(n, 0.0) != 0:
                P0 = P0 + self.p[n] * self.model.direct_phase(n, self.tl, g, self.p)
        Pj = [self.model.direct_phase(n, self.tl, g, self.p) if n in direct else k @ M
              for n, k in zip(self.lin, self.Kj)]
        dg = np.diff(g)
        idx = np.clip(np.searchsorted(self.base, 0.5 * (g[1:] + g[:-1])) - 1, 0, len(split) - 1)

        def panel_rate(P):
            r = np.zeros(len(split))
            np.maximum.at(r, idx, np.abs(np.diff(P, axis=1)).max(axis=0) / dg)
            return r
        self.grid = dict(split=split, g=g, wg=wg, P0=P0, Pj=Pj,
                         R0=panel_rate(P0), Rj=[panel_rate(P) for P in Pj])

    def _phase(self, theta):
        """Total phase, refining the grid when the bound V + Σ|θ_j| R_j requires it."""
        if not hasattr(self, "base"):
            self.base = g_base_edges(self.v, self.acc)
            n = len(self.base) - 1
            self._set_grid(_SPLIT_CACHE.get((self.m, self.w0, self.v, self.acc), (1,) * n))
        width = np.diff(self.base)
        for _ in range(8):
            G = self.grid
            rate = self.v + G["R0"] + sum(abs(th) * R for th, R in zip(theta, G["Rj"]))
            need = np.maximum(rate * width / (2 * np.pi), 1.0)
            need = 2 ** np.ceil(np.log2(np.minimum(need, MAX_SPLIT)))
            if np.all(need <= np.asarray(G["split"])):
                break
            self._set_grid(tuple(int(max(a, b)) for a, b in zip(need, G["split"])))
        G = self.grid
        phi = G["P0"] + self.v * G["g"][None, :]
        for th, P in zip(theta, G["Pj"]):
            phi = phi + th * P
        return phi

    def unit_phase(self):
        """Largest |Φ_j| (in rad) per unit of each θ_j, on the current grid."""
        if not self.ok:
            return [1.0] * len(self.lin)
        if not hasattr(self, "base"):
            self._phase([0.0] * len(self.lin))
        return [float(np.abs(P).max()) or 1.0 for P in self.grid["Pj"]]

    def intensity(self, theta, scale=1.0, jac=False):
        """I(t) for the given linear θ. With jac=True also returns
        ∂I/∂θ_j (columns) and ∂I/∂I0, analytically: I = (1+V²)|A|², A = Σ w e^{iφ}."""
        out = np.ones_like(self.t)
        if not self.ok:
            out[self.live] = np.nan
            return (out, None) if jac else out
        v = self.v
        phi = self._phase(theta)
        wg = self.grid["wg"]
        C, S = np.cos(phi), np.sin(phi)
        Ar, Ai = C @ wg, S @ wg
        f = 1 + v * v
        out[self.live] = (Ar * Ar + Ai * Ai) * f
        if not jac:
            return out * scale
        J = np.zeros((self.t.size, len(theta) + 1))
        for j, P in enumerate(self.grid["Pj"][:len(theta)]):
            dr = -(P * S) @ wg
            di = (P * C) @ wg
            J[self.live, j] = 2 * f * (Ar * dr + Ai * di) * scale
        J[:, -1] = out            # ∂I/∂I0
        return out * scale, J


# ---------------------------------------------------------------------------
# TLLAM — Malacarne et al., Appl. Spectrosc. 65, 99 (2011)
# thermal lens in the sample + thermal lens in the coupled fluid
# ---------------------------------------------------------------------------


def _hh(x):
    """(cosh x - 1)/(x + sinh x), stable form (no overflow)."""
    em = np.expm1(-x)
    return em * em / (-np.expm1(-2 * x) + 2 * x * np.exp(-x))


class TLLAM(Model):
    theta_info = {
        "ths": {"def": "θ sample = ϑ·10⁻⁶,  ϑ = P_e A_e φ / (k_s λ_p)", "unit": "10⁶ K/m",
                "note": "The sample phase is ϑ·l·χ(α,l), with χ (dn/dT, α_T, ν, E, q∥, q⊥) in 10⁻⁶/K. "
                        "It is not the θ in rad of the Shen model: the thickness l and the coefficient χ are left out of it."},
        "thf": {"def": "θ fluid = 2 P_e A_e φ (∂n_f/∂T) / λ_p", "unit": "W/(K·m²)",
                "note": "θ fluid / θ sample = 2·10⁶ k_s (∂n_f/∂T): with k_s known, it gives the fluid's ∂n/∂T."},
    }
    linear = ("ths", "thf")

    def __init__(self):
        super().__init__(
            key="tllam",
            title="TL+TL-Fluid (TLLAM)",
            description="Thermal lens with sample–fluid heat coupling "
                        "(Malacarne et al., Appl. Spectrosc. 65, 99, 2011). "
                        "Equivalent to tllam.c, with the fluid term corrected.",
            presets={
                "Air": {"values": dict(df=2.2e-5, cf=1005, denf=1.18, nf=1.0),
                        "bounds": dict(thf=[-200, 200])},
                "Water": {"values": dict(df=1.45e-7, cf=4182, denf=998, nf=1.33),
                          "bounds": dict(thf=[-5000, 5000])},
            },
            params=[
                Param("d", 5e-7, "D sample", "m²/s", True, 1e-8, 1e-5),
                Param("ths", 4.0, "θ sample", "10⁶ K/m", True, -50, 50),
                Param("thf", 0.0, "θ fluid", "W/(K·m²)", True, -200, 200),
                Param("m", 22.1, "m", ""),
                Param("w0", 57e-6, "w0e", "m"),
                Param("v", 2.07, "V", ""),
                Param("l", 1.85e-3, "L thickness", "m"),
                Param("df", 2.2e-5, "D fluid", "m²/s"),
                Param("cs", 810, "c sample", "J/(kg K)"),
                Param("cf", 1005, "c fluid", "J/(kg K)"),
                Param("dens", 2890, "ρ sample", "kg/m³"),
                Param("denf", 1.18, "ρ fluid", "kg/m³"),
                Param("ns", 1.63, "n sample", ""),
                Param("nf", 1.0, "n fluid", ""),
                Param("E", 72, "E (Young)", ""),
                Param("coefexp", 7.5, "α_T expansion", ""),
                Param("dndts", 5.3, "dn/dT sample", ""),
                Param("nu", 0.23, "ν Poisson", ""),
                Param("qpen", 0.9e-3, "q⊥", ""),
                Param("qpar", 0.09e-3, "q∥", ""),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        u2 = (a * p["w0"]) ** 2
        gauss = np.exp(-u2 / 8)
        K = np.zeros(np.broadcast_shapes(t.shape, a.shape))
        # --- sample (argfasetl) ---
        if p["ths"] != 0:
            l, ns, nf, nu, ce = p["l"], p["ns"], p["nf"], p["nu"], p["coefexp"]
            x = l * a
            hl = _hh(x) / x
            qci2 = ns ** 3 * p["E"] * ce / (4 * (1 - nu))
            B = (p["dndts"] + 4 * (ns - nf) * (1 + nu) * ce * hl
                 + qci2 * (p["qpar"] + 3 * p["qpen"]
                           - 4 * (p["qpar"] * nu + p["qpen"] * (2 + nu)) * hl))
            K = K + p["ths"] * l * gauss * B * (-np.expm1(-p["d"] * t * a * a)) / a
        # --- fluid (argfasear) ---
        if p["thf"] != 0:
            K = K + p["thf"] * self._fluid(t, a, p, gauss)
        return K

    @staticmethod
    def _fluid(t, a, p, gauss):
        d, df = p["d"], p["df"]
        k = d * p["cs"] * p["dens"]
        kf = df * p["cf"] * p["denf"]
        k2, kf2 = k * k, kf * kf
        den = k2 * df - kf2 * d
        a2t = a * a * t
        if den == 0:
            # t0 → ±∞: the terms with t0 go to zero from both sides (continuous limit)
            return gauss / ((k2 - kf2) * a * a) * (
                k * erf(a * np.sqrt(df * t)) - kf * erf(a * np.sqrt(d * t)))
        t0 = (k2 - kf2) * d * df / den
        temp1 = gauss / ((k2 - kf2) * a * a)
        bb1 = k * erf(a * np.sqrt(df * t)) - kf * erf(a * np.sqrt(d * t))

        # Each term is  A·exp(−t0 α² t) + R.  In the erf branch of the original,
        #   −c s e^{−t0α²t} erf(x) = −c s e^{−t0α²t} + c s e^{−dd α² t} erfcx(x),
        # with x² = (dd − t0) α² t. When both terms use this branch,
        # k·s(df) = kf·s(d) exactly and the e^{−t0α²t} parts cancel — which
        # matters for t0 < 0 (e.g. sample in water), where e^{−t0α²t} overflows.
        # Branch t0 > dd: analytic continuation of the erf branch (erf(iy)/i = erfi(y)).
        # The original tllam.c had this term with the wrong sign; checked against
        # the numerical inversion of the Laplace-space solution (tests/validate.py).
        def branch(dd, coef):
            if t0 > dd:
                return 0.0, -(2 * coef / SQRT_PI) * np.sqrt(dd / (t0 - dd)) * \
                    np.exp(-a2t * dd) * dawsn(a * np.sqrt((t0 - dd) * t))
            s = coef * np.sqrt(dd / (dd - t0))
            return -s, s * np.exp(-dd * a2t) * erfcx(a * np.sqrt((dd - t0) * t))

        A3, R3 = branch(d, kf)
        A2, R2 = branch(df, k)
        bb = R2 - R3
        if not (t0 < d and t0 < df):        # only one erf branch: here t0 > 0
            bb = bb + (A2 - A3) * np.exp(-t0 * a2t)
        return temp1 * (bb1 + bb)


# ---------------------------------------------------------------------------
# Shen thermal lens (closed form, mode-mismatched)
#   Shen, Lowe & Snook, Chem. Phys. 165, 385 (1992); Baesso et al., JAP 75, 3732 (1994)
#   I(t) = I(0) [1 − (θ/2) atan(2mV / ([(1+2m)² + V²] t_c/(2t) + 1 + 2m + V²))]²
# ---------------------------------------------------------------------------


class TLShen(Model):
    theta_info = {
        "th": {"def": "θ = −P_e A_e L φ (ds/dT) / (k λ_p)", "unit": "rad",
               "note": "θ of the classical thermal lens: includes the thickness L. Do not confuse it with the θ_TM of the "
                       "thermal mirror (1/m) nor with the θ sample of TLLAM. ds/dT > 0: θ < 0 and the signal rises."},
    }
    linear = ()
    phase_based = False     # no integral phase: the fit uses PSO + LM directly
    proportional = ("th",)  # in a series, θ ∝ P (slope θ/P)

    def __init__(self):
        super().__init__(
            key="tlshen",
            title="Shen TL (classical)",
            description="Closed-form expression of Shen's aberrant model for the thermal "
                        "lens in the mode-mismatched configuration (thin sample, low absorption, small "
                        "θ). t_c = w0e²/(4D). Baesso, Shen & Snook, J. Appl. Phys. "
                        "75, 3732 (1994), eq. 10.",
            params=[
                Param("d", 1.4e-7, "D sample", "m²/s", True, 1e-9, 1e-4),
                Param("th", -0.1, "θ", "rad", True, -3, 3),
                Param("m", 30.0, "m", ""),
                Param("w0", 50e-6, "w0e", "m"),
                Param("v", 2.0, "V", ""),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def intensity(self, t, p, acc=Accuracy(), info=None):
        t = np.asarray(t, float)
        m, v = float(p["m"]), float(p["v"])
        tc = p["w0"] ** 2 / (4 * p["d"])
        with np.errstate(divide="ignore"):
            den = ((1 + 2 * m) ** 2 + v * v) * tc / (2 * np.maximum(t, 0.0)) + 1 + 2 * m + v * v
        arg = np.where(t > 0, 2 * m * v / den, 0.0)
        out = (1 - 0.5 * p["th"] * np.arctan(arg)) ** 2
        if info is not None:
            info.update(n_alpha=0, n_g=0, max_phase=float(abs(p["th"])))
        return out * p.get("scale", 1.0)


# ---------------------------------------------------------------------------
# Thermal lens with photobleaching — Herculano et al., J. Phys. Chem. B 117, 1932 (2013)
#   β(t) = β [(1−ε̄) e^{−K_T t} + ε̄]  (absorption decaying from β to ε̄·β at rate K_T)
#   Φ(g,t) = θ/t_c ∫₀ᵗ [ε̄ + (1−ε̄) e^{−K_T (t−ξ)}] (1+2ξ/t_c)⁻¹ [1 − e^{−2mg/(1+2ξ/t_c)}] dξ
# In Hankel space the same result is  Φ = ∫ K(α,t)[J0(α w0 √(m g)) − 1] dα  with
#   K = −θ e^{−α²w²/8} F(α,t)/α,
#   F = ε̄ (1 − e^{−bt}) + (1−ε̄) b (e^{−bt} − e^{−K_T t})/(K_T − b),   b = D α²
# (convolution of β(t) with the thermal response; K_T → 0 or ε̄ = 1 give Shen's classical TL).
# ---------------------------------------------------------------------------


def _phi1(z):
    """(1 − e^{−z})/z, stable as z → 0."""
    z = np.abs(z)
    small = z < 1e-8
    zs = np.where(small, 1.0, z)
    return np.where(small, 1.0 - 0.5 * z, -np.expm1(-zs) / zs)


class TLFotobranq(Model):
    theta_info = {
        "th": {"def": "θ = −P_e β L φ (dn/dT) / (k λ_p),  β = absorption at t = 0 (before bleaching)",
               "unit": "rad",
               "note": "Same θ as the Shen TL, but with the initial β: the long-time signal corresponds to ε̄·θ. "
                       "dn/dT < 0 (water): θ > 0 and the signal drops. In a series, θ ∝ P_e. "
                       "K_T [1/s] = K_m + K_pb P_e (K_m: mass diffusion; K_pb: photobleaching) "
                       "and ε̄ = ε_products/ε_reactants — the K_T × P_e line appears in the series."},
    }
    linear = ("th",)

    def __init__(self):
        super().__init__(
            key="tlbranq",
            title="TL with photobleaching (Herculano 2013)",
            description="Thermal lens of a liquid whose absorption decays by a photochemical reaction "
                        "(pseudo-first order, rate K_T) down to ε̄·β. Fits θ, K_T and ε̄; fix D with a "
                        "measurement of the pure solvent (K_T = 0). Herculano et al., J. Phys. Chem. B 117, "
                        "1932 (2013), eqs. 3–6.",
            params=[
                Param("d", 1.41e-7, "D sample", "m²/s", False, 1e-9, 1e-4),
                Param("th", 0.05, "θ", "rad", True, -5, 5),
                Param("kt", 2.0, "K_T", "1/s", True, 1e-4, 1e3),
                Param("eb", 0.3, "ε̄", "", True, 0, 5),
                Param("m", 35.1, "m", ""),
                Param("w0", 70e-6, "w0e", "m"),
                Param("v", 7.5, "V", ""),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        b = p["d"] * a * a
        k, eb = p["kt"], p["eb"]
        bt = b * t
        # (e^{−bt} − e^{−kt})/(k − b) = t e^{−min(k,b)t} (1 − e^{−z})/z, z = |k − b| t
        mix = t * np.exp(-np.minimum(k, b) * t) * _phi1((k - b) * t)
        F = eb * (-np.expm1(-bt)) + (1 - eb) * b * mix
        return -p["th"] * np.exp(-(a * p["w0"]) ** 2 / 8) * F / a


# ---------------------------------------------------------------------------
# Thermal lens + population lens — Silva et al., Appl. Phys. Lett. 102, 141910 (2013)
#   Φ = Φ_th + Φ_PL
#   Φ_th = θ_th ∫₀ᵗ dξ/(t_c + 2ξ) [1 − e^{−2mg/(1+2ξ/t_c)}]            (= Shen TL)
#   Φ_PL = θ_PL ∫₀ˡ  u/(1/τ + h u) [1 − e^{−t(1/τ + h u)}] dz,   u = e^{−2mg} e^{−A_e z}
# (excited-state population from a saturable two-level rate equation; h is the
# excitation rate at the centre, h = I0 (σ_abs + σ_se)/E). With u = e^{−2mg − A_e z} the z integral becomes
#   Φ_PL = θ_PL /(A_e h) ∫_{s_L}^{s_0} (1 − e^{−t s})/s ds,  s = 1/τ + h u,
#   ∫ = ln(s0/sL) − [E1(t sL) − E1(t s0)]   (series in t s for t s0 ≤ 1; quadrature in z if A_e l < 0.02).
# As in the paper, the phase is Φ_PL(g) − Φ_PL(0) (a constant phase does not change the intensity).
# ---------------------------------------------------------------------------

from scipy.special import exp1 as _exp1

_GL8 = np.polynomial.legendre.leggauss(8)


def _pl_unit_phase(t, u0, ae, l, h, tau):
    """Φ_PL for θ_PL = 1; t (Nt,1), u0 = e^{−2mg} (1,Ng) or (Ng,)."""
    t = np.asarray(t, float)
    u0 = np.asarray(u0, float)[None, :] if np.ndim(u0) == 1 else u0
    k = 1.0 / tau
    if ae * l < 0.02:
        x, w = _GL8
        out = 0.0
        for xi, wi in zip(x, w):
            z = 0.5 * l * (xi + 1)
            u = u0 * np.exp(-ae * z)
            sk = k + h * u
            out = out + 0.5 * l * wi * u / sk * (-np.expm1(-t * sk))
        return out
    s0 = k + h * u0
    sL = k + h * u0 * np.exp(-ae * l)
    x0 = t * s0
    small = x0 <= 1.0
    # series: Σ (−1)^{j+1} t^j (s0^j − sL^j)/(j·j!)
    a0, aL = np.where(small, x0, 0.0), np.where(small, t * sL, 0.0)
    ser = 0.0
    fact = 1.0
    for j in range(1, 15):
        fact *= j
        ser = ser + (-1) ** (j + 1) * (a0 ** j - aL ** j) / (j * fact)
    xL = np.where(small, 1.0, t * sL)
    xb = np.where(small, 1.0, x0)
    closed = np.log(s0 / sL) - (_exp1(xL) - _exp1(xb))
    return np.where(small, ser, closed) / (ae * h)


class TLPopulacao(Model):
    theta_info = {
        "th": {"def": "θ_th = −P_abs (ds/dT) φ / (k λ_p),  P_abs = P A_e L,  L = (1 − e^{−A_e l})/A_e",
               "unit": "rad",
               "note": "Same θ as Shen. ds/dT < 0 (ZBLAN) with heating: θ_th > 0; "
                       "laser cooling reverses the sign (φ < 0)."},
        "tpl": {"def": "θ_PL = I0 σ_abs N_T Δα_P (2π f_L)² / (n λ_p E)",
                "unit": "1/(m·s)",
                "note": "Strength of the population lens. Δα_P: electronic polarizability difference "
                        "excited–ground state; f_L = (n²+2)/3; E = h c/λ_ex; I0 = 2P/(π w_e²)."},
    }
    linear = ("th", "tpl")
    direct = ("tpl",)

    def __init__(self):
        super().__init__(
            key="tlpop",
            title="TL + population lens (Silva 2013)",
            description="Thermal lens plus excited-state population lens (lifetime τ, excitation "
                        "rate h = I0(σ_abs+σ_se)/E at the centre, absorption A_e, thickness l). "
                        "Silva et al., Appl. Phys. Lett. 102, 141910 (2013), eqs. 1–3. "
                        "Fits θ_th, θ_PL, τ and D; m and V come from the setup, h, A_e and l are known.",
            params=[
                Param("d", 2.5e-7, "D sample", "m²/s", True, 1e-9, 1e-4),
                Param("th", 0.04, "θ_th", "rad", True, -3, 3),
                Param("tpl", 2900.0, "θ_PL", "1/(m·s)", True, -1e7, 1e7),
                Param("tau", 1.7e-3, "τ (lifetime)", "s", True, 1e-5, 1e-1),
                Param("h", 40.0, "h (excitation rate)", "1/s", False, 1e-3, 1e5),
                Param("ae", 240.0, "A_e absorption", "1/m", False, 0.0, 1e5),
                Param("l", 3.79e-3, "l thickness", "m"),
                Param("m", 28.7, "m", ""),
                Param("w0", 144e-6, "w0e", "m"),
                Param("v", 15.0, "V", ""),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        # thermal part (Shen): K = −θ_th e^{−α²w²/8} (1 − e^{−Dα²t})/α
        return -p["th"] * np.exp(-(a * p["w0"]) ** 2 / 8) * (-np.expm1(-p["d"] * a * a * t)) / a

    def direct_phase(self, name, t, g, p):
        g = np.asarray(g, float)
        u0 = np.exp(-2 * p["m"] * g)
        both = np.concatenate([[1.0], u0])          # u = 1 at g = 0 (centre reference)
        ph = _pl_unit_phase(t, both, p["ae"], p["l"], p["h"], p["tau"])
        return ph[:, 1:] - ph[:, :1]


# ---------------------------------------------------------------------------
# Thermal mirror: LAM, BAM (Beer's law) and LAM with finite thickness
#   Sato et al., J. Appl. Phys. 104, 053520 (2008)  — LAM, BAM, HAM
#   Astrath et al., J. Opt. Soc. Am. B 28, 1735 (2011) — finite thickness
#
# With θ_TM = −P_e A_e α_T (1+ν) φ/(λ_p k) [1/m], the phase is
#   Φ(g,t) = 4 θ_TM ∫ e^{−α²w²/8} F(α,t) [J0(α w √(mg)) − 1] dα,
#   F_LAM  = (1 − e^{−Dα²t})/α²                       (Sato eq. 24/32)
#   F_fin  = h(lα) (1 − e^{−Dα²t})/α²,  h(x) = (cosh x − 1)/(x + sinh x)
#                                                      (Astrath 2011, eq. 10)
#   F_BAM  = α D Ĝ(α,t),  Ĝ = ∫_0^∞ T̂(α,z,t) e^{−αz} dz  (source e^{−A_e z})
# Ĝ comes from dĜ/dt = [A_e erfc(αs) − α e^{−α²s²} erfcx(A_e s)]/(A_e² − α²),
# s = √(Dt), integrated in closed form (or by a series when s·max(α, A_e) ≤ 0.75).
# Limits: A_e → 0 recovers LAM; A_e → ∞ recovers HAM (θ_HAM = θ_TM/A_e).
# Switch-off at t_off by superposition: F(t) − F(t − t_off).
# ---------------------------------------------------------------------------

from math import factorial as _factorial
from scipy.special import gamma as _gamma_fn

_XERFC_C = np.array([(-1) ** n / (_factorial(n) * (2 * n + 1) * (2 * n + 3)) for n in range(14)])


def _xerfc_int(X):
    """∫_0^X x erfc(x) dx (Taylor series for small X, no cancellation)."""
    out = np.empty_like(X)
    sm = X < 0.3
    xs = X[sm]
    out[sm] = xs ** 2 / 2 - (2 / SQRT_PI) * xs ** 3 * np.polyval(_XERFC_C[::-1], xs ** 2)
    xb = X[~sm]
    out[~sm] = xb * xb * erfc(xb) / 2 + erf(xb) / 4 - xb * np.exp(-xb * xb) / (2 * SQRT_PI)
    return out


def _bam_closed(a, t, ae, d):
    s = np.sqrt(d * t)
    x, y = a * s, ae * s
    i1 = 2 * t * _xerfc_int(x) / (x * x)                     # ∫_0^t erfc(α√(Dτ)) dτ
    c = (ae * ae - a * a) * d
    i2 = (np.exp(-x * x) * erfcx(y) - 1) / c + ae * erf(x) / (c * a)
    return (ae * i1 - a * i2) / (ae * ae - a * a)


_BAM_NS = 30
_BAM_GH = _gamma_fn(np.arange(_BAM_NS) / 2 + 1)
_BAM_KF = np.array([1 / _factorial(k) for k in range(_BAM_NS // 2)])


def _bam_series(a, t, ae, d):
    """Ĝ = t/(A_e+α) Σ_n p_n ŝ^n/(n/2+1), in variables scaled by σ = max(α, A_e).

    Comes from erfcx(z) = Σ (−z)^n/Γ(n/2+1): the factor (A_e² − α²) cancels term by
    term, so there is no division by small numbers near α = A_e. The coefficients p_n
    depend only on α: they are computed once per distinct α."""
    au, inv = np.unique(a, return_inverse=True)
    sig = np.maximum(au, ae)
    ah, Ah = au / sig, ae / sig
    n = _BAM_NS
    S = np.empty((n, au.size))
    S[0] = 1.0
    ahk = np.ones(au.size)
    for k in range(1, n):
        ahk *= ah
        S[k] = Ah * S[k - 1] + ahk
    A = np.zeros((n, au.size))
    A[0] = 1.0
    c = -Ah * ah
    for k in range(2, n):
        A[k] = ((-1) ** k / _BAM_GH[k]) * c * S[k - 2]
    E = np.zeros((n, au.size))
    E[0] = 1.0
    q = -(ah * ah)
    for j in range(1, n // 2):
        E[2 * j] = E[2 * j - 2] * q / j
    P = np.zeros((n, au.size))
    for k in range(0, n, 2):                 # E has only even terms
        P[k:] += E[k] * A[:n - k]
    P /= (np.arange(n) / 2 + 1)[:, None]
    P = P[:, inv]
    sh = sig[inv] * np.sqrt(d * t)
    acc = np.zeros(a.size)
    for k in range(n - 1, -1, -1):
        acc = acc * sh + P[k]
    return t / (ae + a) * acc


def bam_g(a, t, ae, d, rel=1e-4, zmax=0.75):
    """Ĝ(α,t) of the BAM model (source e^{−A_e z}, free and adiabatic surface).

    Checked against direct quadrature of dĜ/dt (relative error < 3e-7 for
    10⁻¹ < A_e < 10⁸ m⁻¹, 10⁻⁶ < t < 10 s, including α ≈ A_e)."""
    shape = np.broadcast_shapes(np.shape(a), np.shape(t))
    a = np.broadcast_to(np.asarray(a, float), shape).ravel()
    t = np.broadcast_to(np.asarray(t, float), shape).ravel()
    out = np.zeros(a.size)
    pos = t > 0
    z = np.maximum(a, ae) * np.sqrt(d * np.where(pos, t, 0.0))
    ser = pos & (z <= zmax)
    if np.any(ser):
        out[ser] = _bam_series(a[ser], t[ser], ae, d)
    cl = pos & ~ser
    if np.any(cl):
        ac, tcl = a[cl], t[cl]
        with np.errstate(all="ignore"):
            g = _bam_closed(ac, tcl, ae, d)
        near = np.abs(ac / ae - 1) < rel
        if np.any(near):   # removable singularity at α = A_e: linear interpolation
            tn = tcl[near]
            lo = _bam_closed(np.full(tn.size, ae * (1 - rel)), tn, ae, d)
            hi = _bam_closed(np.full(tn.size, ae * (1 + rel)), tn, ae, d)
            g[near] = 0.5 * (lo + hi) + (ac[near] - ae) / (2 * rel * ae) * (hi - lo)
        out[cl] = g
    return out.reshape(shape)


class _TMBase(Model):
    """Thermal mirror (phase 4 θ_TM ∫ e^{−α²w²/8} F(α,t) [J0 − 1] dα) with on/off."""
    linear = ("th",)
    theta_info = {
        "th": {"def": "θ_TM = −P_e A_e α_T (1+ν) φ / (λ_p k)", "unit": "1/m",
               "note": "Same θ_TM in the LAM, BAM and finite-sample models (it is the θ_TM of the papers, in m⁻¹). "
                       "The HAM model uses θ_TM/A_e, dimensionless. Expanding sample (α_T > 0): θ_TM < 0."},
    }

    def _F(self, t, a, p):
        raise NotImplementedError

    def kernel(self, t, a, p, acc):
        tt = np.broadcast_to(t, np.broadcast_shapes(np.shape(t), np.shape(a)))
        F = self._F(tt, a, p)
        of = p.get("of", 0.0)
        if of and of > 0 and np.any(tt > of):
            F = F - self._F(np.maximum(tt - of, 0.0), a, p)
        return 4 * p["th"] * np.exp(-(a * p["w0"]) ** 2 / 8) * F


_TM_COMMON = lambda th_lo, th_hi: [
    Param("d", 6e-7, "D sample", "m²/s", True, 1e-8, 1e-4),
    Param("th", -300.0, "θ TM", "1/m", True, th_lo, th_hi),
    Param("m", 20.0, "m", ""),
    Param("w0", 44e-6, "w0e", "m"),
    Param("v", 5.2, "V", ""),
]


class TMLAM(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmlam",
            title="TM LAM (low absorption)",
            description="Thermal mirror for a semi-infinite sample with a source uniform "
                        "in z (A_e·thickness ≪ 1). LAM = Low Absorption Model. Sato et al., J. Appl. Phys. 104, "
                        "053520 (2008), eq. 32. θ_TM = −P_e A_e α_T(1+ν)φ/(λ_p k), in 1/m.",
            params=_TM_COMMON(-1e6, 1e6) + [
                Param("of", 0.0, "t_off (0 = no switch-off)", "s"),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        return -np.expm1(-p["d"] * a * a * t) / (a * a)


class TMFinito(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmfinito",
            title="TM LAM, finite sample",
            description="Low-absorption thermal mirror in a sample of thickness l, "
                        "with both faces free. Astrath et al., J. Opt. Soc. Am. B 28, "
                        "1735 (2011), eq. 10. Tends to the semi-infinite LAM for large l.",
            params=_TM_COMMON(-1e6, 1e6) + [
                Param("l", 1e-3, "l thickness", "m"),
                Param("of", 0.0, "t_off (0 = no switch-off)", "s"),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        return _hh(p["l"] * a) * -np.expm1(-p["d"] * a * a * t) / (a * a)


class TMBAM(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmbam",
            title="TM BAM (Beer's law)",
            description="Thermal mirror with absorption following Beer's law, valid for "
                        "any A_e (semi-infinite sample); BAM = Beer Absorption Model. Sato et al., J. Appl. Phys. "
                        "104, 053520 (2008), eq. 33–34. Tends to LAM (A_e → 0) and to HAM "
                        "(A_e → ∞, with θ_HAM = θ_TM/A_e).",
            params=_TM_COMMON(-1e9, 1e9) + [
                Param("ae", 1000.0, "A_e absorption", "1/m", False, 1.0, 1e8),
                Param("of", 0.0, "t_off (0 = no switch-off)", "s"),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        d = p["d"]
        return a * d * bam_g(a, t, float(p["ae"]), d)


# ---------------------------------------------------------------------------
# TM+TL-Fluid, HAM on/off (approximation, no feedback from the fluid) — equivalent to tmhamtl-onoff.c
# ---------------------------------------------------------------------------


class TMHAMTL(Model):
    theta_info = {
        "th": {"def": "θ_HAM = θ_TM / A_e = −P_e α_T (1+ν) φ / (λ_p k)", "unit": "dimensionless",
               "note": "Does not contain A_e (opaque sample). To compare with LAM, BAM and finite sample, "
                       "which give θ_TM in 1/m: θ_TM = θ_HAM·A_e. Expanding sample (α_T > 0): θ < 0. "
                       "With incident P_e and reflectivity R, the θ_TM of Flizikowski et al. (2020) is (1−R)·θ_HAM."},
    }
    linear = ("th",)     # the air lens is also proportional to θ (pref ∝ th)

    def __init__(self):
        super().__init__(
            key="tmhamtl",
            title="TM+TL-Fluid, HAM on/off (approximation: no feedback from the fluid)",
            description="HAM thermal mirror (High Absorption Model) with the excitation switched off at "
                        "t = t_off, plus the thermal lens of the fluid (air), active if dn/dT ≠ 0. "
                        "Equivalent to tmhamtl-onoff.c: approximation in which the sample surface is "
                        "free; the fluid receives heat and forms the TL, but does not cool the sample back. "
                        "For the full coupling (Flizikowski et al., Opt. Express 28, 7116, 2020) "
                        "use TM+TL-Fluid, HAM (full coupling): in the three metals of the paper the "
                        "difference is 1e-4 to 1e-3 of the signal with the laser on and reaches ~1% right "
                        "after switch-off (stainless steel in air).",
            params=[
                Param("d", 3.46e-6, "D sample", "m²/s", True, 1e-8, 1e-4),
                Param("th", -0.57, "θ HAM", "dimensionless", True, -5, 5),
                Param("m", 6.01, "m", ""),
                Param("w0", 646e-6, "w0e", "m"),
                Param("v", 52.6, "V", ""),
                Param("beta", 1.73e-5, "β expansion", "1/K"),
                Param("nu", 0.29, "ν Poisson", ""),
                Param("dndt", 0.0, "dn/dT air", "1/K"),
                Param("df", 2.2e-5, "D air", "m²/s"),
                Param("of", 0.899, "t_off", "s"),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        d, th, w0, of = p["d"], p["th"], p["w0"], p["of"]
        u2 = (a * w0) ** 2
        gauss = np.exp(-u2 / 8)
        a2 = a * a
        tt = np.broadcast_to(t, np.broadcast_shapes(t.shape, a.shape))
        on = tt <= of
        # argfasetm (t <= of)
        st = np.sqrt(d * t)
        k_on = 2 * th * gauss * (-2 * np.sqrt(d * t / np.pi) * np.exp(-d * t * a2)
                                 + erf(st * a) / a + erfc(st * a) * 2 * d * t * a)
        K = np.where(on, k_on, 0.0)
        if np.any(~on):
            to = np.maximum(t - of, 0.0)
            sto = np.sqrt(d * to)
            k_off = (2 * np.sqrt(d * to / np.pi) * np.exp(-a2 * (d * to + w0 * w0 / 8))
                     - 2 * np.sqrt(d * t / np.pi) * np.exp(-a2 * (d * t + w0 * w0 / 8))
                     + gauss / a * ((2 * d * t * a2 - 1) * erfc(st * a)
                                    - (2 * d * to * a2 - 1) * erfc(sto * a)))
            K = np.where(on, K, 2 * th * k_off)
        if p["dndt"] != 0:
            K = K + self._air(t, a, p, gauss, acc.tau_nodes)
        return K

    @staticmethod
    def _air(t, a, p, gauss, n_tau=80, block=16):
        """argfasear/argtau: τ integral done with τ = t sin²φ and tanh-sinh.

        Vectorized over (time × α × φ), in blocks of `block` times to limit
        memory (~block·Nα·Nφ numbers per block)."""
        d, df, of = p["d"], p["df"], p["of"]
        pref = -p["th"] * 2 * p["dndt"] / (p["beta"] * (1 + p["nu"])) * np.sqrt(df / np.pi)
        x, w = tanh_sinh(n_tau)
        tv = t[:, 0]
        av2 = (a[0] ** 2)[None, :, None]                      # (1, Nα, 1)
        av = a[0][None, :, None]
        J = np.empty((tv.size, a.shape[1]))
        for i0 in range(0, tv.size, block):
            ti = tv[i0:i0 + block]
            lo = np.where(ti <= of, 0.0, np.arcsin(np.sqrt(np.clip((ti - of) / ti, 0, 1))))
            span = 0.5 * np.pi - lo
            phi = lo[:, None] + span[:, None] * x[None, :]   # (nb, Nφ)
            wp = span[:, None] * w[None, :]
            s2 = np.sin(phi) ** 2
            c = np.cos(phi)
            tt = ti[:, None, None]
            # 2 sqrt(t) ∫ cosφ e^{-α² df t sin²φ} erf(α sqrt(d t) cosφ) dφ
            f = c[:, None, :] * np.exp(-av2 * df * tt * s2[:, None, :]) * \
                erf(av * np.sqrt(d * tt) * c[:, None, :])
            J[i0:i0 + block] = 2 * np.sqrt(ti)[:, None] * np.einsum("tap,tp->ta", f, wp)
        return pref * gauss * J


# ---------------------------------------------------------------------------
# Thermal mirror (opaque sample) with thermal lens in the coupled fluid — with full
# sample–fluid heat coupling.
#   Lukasievicz et al., Appl. Spectrosc. 66, 1461 (2012), eqs. 4–9, 18–22
#
# Surface source; sample at z>0 and fluid at z<0. In Hankel–Laplace (q_i² = α² + s/D_i)
#   T(α, z=0, s) ∝ 1 / [s (k_s q_s + k_f q_f)],
#   Ĝ_s = ∫_0^∞ T e^{−αz} dz ∝ 1/[s (k_s q_s + k_f q_f)(q_s + α)]   (mirror)
#   Ĝ_f = ∫_{−∞}^0 T dz     ∝ 1/[s q_f (k_s q_s + k_f q_f)]         (lens in the fluid)
# With ŝ = s/(D_s α²), τ = D_s α² t, κ = k_f/k_s and r = D_s/D_f, both depend only on τ:
#   g_s(τ) = L⁻¹{1/[ŝ (√(1+ŝ) + κ√(1+rŝ)) (√(1+ŝ) + 1)]}
#   g_f(τ) = L⁻¹{1/[ŝ √(1+rŝ) (√(1+ŝ) + κ√(1+rŝ))]}
# tabulated in log τ by numerical Laplace inversion (Talbot) and interpolated.
# Phase  Φ = ∫ K [J0(α w0 √(mg)) − 1] dα,  K = e^{−α²w²/8} (4 h_TM g_s + (h_f/k_s) g_f)/α.
# Without fluid (κ = 0) K_TM is that of HAM (θ_HAM = h_TM).
# ---------------------------------------------------------------------------

from scipy.interpolate import CubicSpline as _CubicSpline

_TAU_LO, _TAU_HI, _TAU_N = -12.0, 9.5, 2600
_TM_TABLES: dict = {}


def _talbot(F, tau, M=22):
    """Inverse Laplace transform (Abate–Valkó, fixed Talbot) of F(ŝ) at several τ; ~1e-12."""
    r = 2 * M / (5 * tau)
    th = np.arange(1, M) * np.pi / M
    cot = 1 / np.tan(th)
    s = r[None, :] * th[:, None] * (cot[:, None] + 1j)
    sig = th + (th * cot - 1) * cot
    f0 = 0.5 * np.real(F(r.astype(complex))) * np.exp(r * tau)
    terms = np.real(np.exp(tau[None, :] * s) * F(s) * (1 + 1j * sig[:, None]))
    return (r / M) * (f0 + terms.sum(0))


def _tm_table(kind, kappa, rr):
    """Spline of g_s ('s') or g_f ('f') in log10 τ, for given (κ, r)."""
    key = (kind, float(kappa), float(rr))
    hit = _TM_TABLES.get(key)
    if hit is not None:
        return hit
    lt = np.linspace(_TAU_LO, _TAU_HI, _TAU_N)
    tau = 10.0 ** lt
    if kind == "s":
        F = lambda x: 1 / (x * (np.sqrt(1 + x) + kappa * np.sqrt(1 + rr * x)) * (np.sqrt(1 + x) + 1))   # noqa: E731
    else:
        F = lambda x: 1 / (x * np.sqrt(1 + rr * x) * (np.sqrt(1 + x) + kappa * np.sqrt(1 + rr * x)))    # noqa: E731
    g = _talbot(F, tau)
    if len(_TM_TABLES) > 48:
        _TM_TABLES.pop(next(iter(_TM_TABLES)))
    _TM_TABLES[key] = _CubicSpline(lt, g)
    return _TM_TABLES[key]


def _tm_g(kind, tau, kappa, rr):
    """g_s or g_f at any τ (asymptotes outside the table)."""
    tau = np.asarray(tau, float)
    out = np.empty_like(tau)
    lo, hi = tau < 10.0 ** _TAU_LO, tau > 10.0 ** _TAU_HI
    mid = ~(lo | hi)
    sq = np.sqrt(rr)
    if kind == "s":
        out[lo] = tau[lo] / (1 + kappa * sq)
        out[hi] = 0.5 / (1 + kappa)
    else:
        out[lo] = tau[lo] / (sq * (1 + kappa * sq))
        out[hi] = 1.0 / (1 + kappa)
    out[mid] = _tm_table(kind, kappa, rr)(np.log10(tau[mid]))
    return out


class TMFluido(Model):
    theta_info = {
        "th": {"def": "h_TM = θ_HAM = −P_e α_T (1+ν) φ / (λ_p k_s)", "unit": "dimensionless",
               "note": "Same θ as TM-HAM (opaque sample): does not contain A_e. Expanding sample: h_TM < 0. "
                       "With incident P_e and reflectivity R, the θ_TM of Flizikowski et al. (2020) is (1−R)·h_TM."},
        "thf": {"def": "h_f = 2 P_e φ (∂n/∂T)_f / λ_p", "unit": "W/(K·m)",
                "note": "Thermal lens of the fluid (two passes: before and after the reflection). "
                        "With k_s known, h_f/P_e gives (∂n/∂T)_f: h_f = 2 P_e φ (∂n/∂T)_f/λ_p. "
                        "P_e is the absorbed power; with incident P_e use P_e(1−R)."},
    }
    linear = ("th", "thf")

    def __init__(self):
        super().__init__(
            key="tmfluido",
            title="TM+TL-Fluid, HAM (full coupling, Flizikowski 2020)",
            description="Thermal mirror of an opaque sample and thermal lens of the non-absorbing fluid "
                        "in contact, with sample–fluid heat exchange (without the restriction k_f√D ≪ k√D_f). "
                        "Lukasievicz et al., Appl. Spectrosc. 66, 1461 (2012), eqs. 20–22. "
                        "k_s = D_s ρ_s c_s and k_f = D_f ρ_f c_f: fix ρ and c and fit the diffusivities. "
                        "In air, the fluid barely changes the sample; in water, it does. "
                        "It is the full model of Flizikowski et al., Opt. Express 28, 7116 (2020), eqs. 3–10; "
                        "TM+TL-Fluid, HAM on/off (approximation) neglects the feedback of the fluid on the sample.",
            presets={
                "Air": {"values": dict(df=2.19e-5, cf=1005, denf=1.18),
                        "bounds": dict(thf=[-1e3, 1e3], df=[1e-6, 1e-3])},
                "Water": {"values": dict(df=1.45e-7, cf=4182, denf=998),
                          "bounds": dict(thf=[-1e5, 1e5], df=[1e-8, 1e-5])},
            },
            params=[
                Param("d", 6.94e-6, "D sample", "m²/s", True, 1e-8, 1e-4),
                Param("th", -0.5, "h_TM (θ HAM)", "dimensionless", True, -50, 50),
                Param("thf", -100.0, "h_f", "W/(K·m)", True, -1e5, 1e5),
                Param("df", 1.45e-7, "D fluid", "m²/s", False, 1e-8, 1e-3),
                Param("m", 32.8, "m", ""),
                Param("w0", 74e-6, "w0e", "m"),
                Param("v", 6.9, "V", ""),
                Param("cs", 460, "c sample", "J/(kg K)"),
                Param("dens", 7800, "ρ sample", "kg/m³"),
                Param("cf", 4182, "c fluid", "J/(kg K)"),
                Param("denf", 998, "ρ fluid", "kg/m³"),
                Param("of", 0.0, "t_off (0 = no switch-off)", "s"),
                Param("scale", 1.0, "I0 (scale)", "", False, 0.8, 1.2),
            ],
        )

    def _K(self, t, a, p):
        d, df = p["d"], p["df"]
        ks = d * p["cs"] * p["dens"]
        kap = df * p["cf"] * p["denf"] / ks
        rr = d / df
        tau = d * a * a * t
        K = 0.0
        if p["th"] != 0:
            K = K + 4 * p["th"] * _tm_g("s", tau, kap, rr)
        if p["thf"] != 0:
            K = K + (p["thf"] / ks) * _tm_g("f", tau, kap, rr)
        return K / a

    def kernel(self, t, a, p, acc):
        tt = np.broadcast_to(t, np.broadcast_shapes(np.shape(t), np.shape(a)))
        aa = np.broadcast_to(a, tt.shape)
        K = self._K(tt, aa, p)
        of = p.get("of", 0.0)
        if of and of > 0 and np.any(tt > of):
            K = K - np.where(tt > of, self._K(np.maximum(tt - of, 0.0), aa, p), 0.0)
        return np.exp(-(aa * p["w0"]) ** 2 / 8) * K



# Order in the interface: first the TL models, then the TM ones (the first is the default on opening).
MODELS = {m.key: m for m in (TLLAM(), TLShen(), TLFotobranq(), TLPopulacao(),
                             TMLAM(), TMFinito(), TMBAM(), TMHAMTL(), TMFluido())}

# Aliases of preset names (example files and older versions used the Portuguese names)
PRESET_ALIASES = {"ar": "Air", "agua": "Water", "água": "Water"}
