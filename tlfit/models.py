"""
Modelos de lente térmica / espelho térmico resolvidos no tempo.

Reimplementação vetorizada dos modelos de TMTL (tllam.c e tmhamtl-onoff.c).
As fórmulas dos integrandos são exatamente as do código C original; o que muda
é a forma de integrar:

  I(t) = (1+v^2) | ∫_0^Gmax e^{-g} exp(i[v g + Φ(g,t)]) dg |^2
  Φ(g,t) = ∫ K(α,t) [J0(α w0 sqrt(m g)) - 1] dα

O integrando de Φ é separável: K(α,t) depende do tempo e dos parâmetros
térmicos; o fator de Bessel depende só de (α, g, m, w0). Usando quadraturas
fixas (Gauss-Legendre por painéis) em α e em g, Φ vira um produto de matrizes

  Φ[t, g] = K[t, α] @ M[α, g],   M = w_α (J0 - 1)

e M é calculada uma única vez por geometria (m, w0). A parte real e a
imaginária da amplitude saem da mesma soma complexa. No código original
cada Φ era uma integral adaptativa recalculada dentro de outra integral
adaptativa (e duas vezes: uma para cos, outra para sin).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy.special import erf, erfc, erfcx, dawsn, j0, expit

SQRT_PI = np.sqrt(np.pi)

# ---------------------------------------------------------------------------
# Quadraturas
# ---------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _gl(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def composite_gl(edges, n):
    """Gauss-Legendre de n pontos em cada painel [edges[i], edges[i+1]]."""
    x, w = _gl(n)
    edges = np.asarray(edges, float)
    a, b = edges[:-1, None], edges[1:, None]
    nodes = 0.5 * (b - a) * x[None, :] + 0.5 * (b + a)
    weights = 0.5 * (b - a) * w[None, :]
    return nodes.ravel(), weights.ravel()


@dataclass(frozen=True)
class Accuracy:
    """Controle de resolução das quadraturas (fator 1 = padrão)."""
    u_max: float = 14.0        # corte em u = α w0 (e^{-u²/8} ≈ 2e-11)
    g_max: float = 30.0        # mesmo limite do código original
    nodes_per_panel: int = 12
    u_panels_per_period: float = 1.0
    g_panels_per_period: float = 1.0
    u_grading_levels: int = 14  # painéis geométricos perto de u = 0
    tau_nodes: int = 80         # tanh-sinh para a integral em τ (ar, dn/dT)

    def scaled(self, f: float) -> "Accuracy":
        return Accuracy(self.u_max, self.g_max, self.nodes_per_panel,
                        self.u_panels_per_period * f, self.g_panels_per_period * f,
                        self.u_grading_levels, int(self.tau_nodes * f))

    def search(self) -> "Accuracy":
        """Versão mais leve, para a fase de busca do ajuste.

        g até 18 (o termo desprezado é ~2e^{-18} ≈ 3e-8 relativo) e malhas 0,7×:
        erro em I(t) ≲ 3e-7, muito abaixo do ruído típico (~1e-4), com 2–3× menos
        trabalho. O refinamento final usa sempre a precisão normal."""
        return Accuracy(self.u_max, min(self.g_max, 18.0), self.nodes_per_panel,
                        self.u_panels_per_period * 0.7, self.g_panels_per_period * 0.7,
                        self.u_grading_levels, max(40, int(self.tau_nodes * 0.7)))


def u_grid(m: float, w0: float, acc: Accuracy):
    u_top = min(acc.u_max, 2e5 * w0)  # o original integra α em [0, 2e5]
    # período de J0(u sqrt(m g)) em u no pior caso (g = g_max)
    period = 2 * np.pi / np.sqrt(max(m, 1e-12) * acc.g_max)
    width = period / acc.u_panels_per_period
    first = min(width, u_top)
    graded = first * 2.0 ** -np.arange(acc.u_grading_levels, 0, -1)
    uniform = np.linspace(first, u_top, max(1, int(np.ceil((u_top - first) / width))) + 1)
    edges = np.concatenate([[0.0], graded, uniform])
    edges = np.unique(edges)
    return composite_gl(edges, acc.nodes_per_panel)


def g_base_edges(v: float, acc: Accuracy):
    """Painéis-base em g: largura = um período de e^{i v g} (ou 1, se v < 2π)."""
    rate = max(abs(v), 2 * np.pi)
    width = 2 * np.pi / rate / acc.g_panels_per_period
    n = int(np.ceil(acc.g_max / width))
    return np.linspace(0, acc.g_max, n + 1)


def g_grid(v: float, acc: Accuracy, split: tuple):
    """Painéis-base subdivididos (split[i] = 1, 2, 4, ...) onde a fase gira rápido."""
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


MAX_SPLIT = 64  # teto de subpainéis por painel-base (fase > ~400 rad: fora do físico)
_SPLIT_CACHE: dict = {}   # (m, w0, v, acc) → subdivisão usada na última avaliação


def needed_split(phi, g, base_edges, n_panels):
    """Subpainéis necessários em cada painel-base para que a fase total
    (v g + Φ) gire no máximo ~2π por painel. Potências de 2 (reaproveita cache)."""
    with np.errstate(invalid="ignore"):
        rate = np.abs(np.diff(phi, axis=1)).max(axis=0) / np.diff(g)   # entre nós vizinhos
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
    """Escolhe a malha em g e devolve (φ, g, pesos).

    phase_fn(M, g) devolve a fase total (Nt, Ng) na malha dada. Começa pela
    subdivisão usada na avaliação anterior com a mesma geometria (quase sempre
    já serve, e então basta um produto K·M) e refina onde for preciso."""
    key = (m, w0, v, acc)
    base = g_base_edges(v, acc)
    n = len(base) - 1
    split = _SPLIT_CACHE.get(key, (1,) * n)
    for _ in range(8):
        M, g, wg = _bessel_matrix(m, w0, acc, v, split)
        phi = phase_fn(M, g)
        need = needed_split(phi, g, base, n)
        if all(a <= b for a, b in zip(need, split)):
            _SPLIT_CACHE[key] = need      # pode encolher na próxima avaliação
            return phi, g, wg
        split = tuple(max(a, b) for a, b in zip(need, split))
    return phi, g, wg


def refine_split(phi, g, v, base_edges, split, nodes_per_panel):
    """Compatibilidade: subdivisão necessária, nunca menor que a atual."""
    need = needed_split(phi, g, base_edges, len(split))
    return tuple(max(a, b) for a, b in zip(need, split))


def tanh_sinh(n: int):
    """Nós/pesos tanh-sinh em [0, 1] (bons para camadas-limite nas pontas)."""
    h = 6.0 / n
    k = np.arange(-n // 2, n // 2 + 1) * h
    s = 0.5 * np.pi * np.sinh(k)
    x = expit(2 * s)
    w = 0.5 * h * 0.5 * np.pi * np.cosh(k) / np.cosh(s) ** 2
    keep = (x > 0) & (x < 1) & (w > 1e-300)
    return x[keep], w[keep]


# ---------------------------------------------------------------------------
# Parâmetros
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
    # predefinições: nome → {"values": {param: valor}, "bounds": {param: [mín, máx]}}
    presets: dict = field(default_factory=dict)

    def defaults(self) -> dict[str, float]:
        return {p.name: p.value for p in self.params}

    # parâmetros que entram linearmente na fase: Φ = Σ θ_j Φ_j (ver PhaseCache).
    # Atributo de classe (sem anotação, para não virar campo da dataclass).
    linear = ()

    # implementado pelas subclasses
    def kernel(self, t, alpha, p, acc):
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
            # parâmetros fora do domínio do modelo: devolve NaN (o ajuste descarta)
            out[live] = np.nan
            if info is not None:
                info.update(n_alpha=len(u), n_g=0, max_phase=float("nan"))
            return out * scale
        phi, g, wg = solve_phase(lambda M, g: K @ M + v * g[None, :], m, w0, v, acc)
        amp = np.cos(phi) @ wg + 1j * (np.sin(phi) @ wg)
        out[live] = (amp.real ** 2 + amp.imag ** 2) * (1 + v * v)
        if info is not None:
            info.update(n_alpha=len(u), n_g=len(g), max_phase=float(np.abs(phi - v * g).max()))
        return out * scale


class PhaseCache:
    """Fase decomposta nos parâmetros lineares, para parâmetros não lineares fixos.

    Φ(t, g) = Φ_0 + Σ_j θ_j Φ_j, com Φ_j = K_j · M e K_j o núcleo com θ_j = 1 e
    os demais θ lineares = 0. Com isso, variar os θ (e a escala I0) custa só o
    cos/sin e a soma em g, 10–15× menos que uma avaliação completa do modelo."""

    def __init__(self, model: Model, t, p: dict, lin_names, acc: Accuracy):
        self.t = np.asarray(t, float)
        self.live = self.t > 1e-9
        self.lin = list(lin_names)
        self.p = dict(p)
        self.acc = acc
        self.m, self.w0, self.v = float(p["m"]), float(p["w0"]), float(p["v"])
        u, _ = _alpha_grid(self.m, self.w0, acc)
        tl = self.t[self.live][:, None]
        a = (u / self.w0)[None, :]
        zero = {n: 0.0 for n in model.linear}
        with np.errstate(all="ignore"):
            # parte fixa: θ lineares que NÃO estão sendo variados, com seus valores
            p0 = {**p, **{n: 0.0 for n in self.lin}}
            self.K0 = model.kernel(tl, a, p0, acc) if any(p0[n] != 0 for n in model.linear) \
                else np.zeros((tl.shape[0], a.shape[1]))
            self.Kj = [model.kernel(tl, a, {**p, **zero, n: 1.0}, acc) for n in self.lin]
        self.ok = bool(np.all(np.isfinite(self.K0)) and all(np.all(np.isfinite(k)) for k in self.Kj))
        self._grid = None   # (id da malha, Φ_0, [Φ_j], g, wg)

    def _set_grid(self, split):
        """Φ_0, Φ_j e as taxas máximas de variação de cada Φ_j por painel-base."""
        M, g, wg = _bessel_matrix(self.m, self.w0, self.acc, self.v, split)
        P0 = self.K0 @ M
        Pj = [k @ M for k in self.Kj]
        dg = np.diff(g)
        idx = np.clip(np.searchsorted(self.base, 0.5 * (g[1:] + g[:-1])) - 1, 0, len(split) - 1)

        def panel_rate(P):
            r = np.zeros(len(split))
            np.maximum.at(r, idx, np.abs(np.diff(P, axis=1)).max(axis=0) / dg)
            return r
        self.grid = dict(split=split, g=g, wg=wg, P0=P0, Pj=Pj,
                         R0=panel_rate(P0), Rj=[panel_rate(P) for P in Pj])

    def _phase(self, theta):
        """Fase total, refinando a malha quando o limite V + Σ|θ_j| R_j exigir."""
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

    def intensity(self, theta, scale=1.0, jac=False):
        """I(t) para os θ lineares dados. Com jac=True devolve também
        ∂I/∂θ_j (colunas) e ∂I/∂I0, analíticos: I = (1+V²)|A|², A = Σ w e^{iφ}."""
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
# lente térmica na amostra + lente térmica no fluido acoplado
# ---------------------------------------------------------------------------


def _hh(x):
    """(cosh x - 1)/(x + sinh x), forma estável (sem overflow)."""
    em = np.expm1(-x)
    return em * em / (-np.expm1(-2 * x) + 2 * x * np.exp(-x))


class TLLAM(Model):
    linear = ("ths", "thf")

    def __init__(self):
        super().__init__(
            key="tllam",
            title="TL amostra + fluido (TLLAM)",
            description="Lente térmica com acoplamento térmico amostra–fluido "
                        "(Malacarne et al., Appl. Spectrosc. 65, 99, 2011). "
                        "Equivalente a tllam.c, com o termo do fluido corrigido.",
            presets={
                "Ar": {"values": dict(df=2.2e-5, cf=1005, denf=1.18, nf=1.0),
                       "bounds": dict(thf=[-200, 200])},
                "Água": {"values": dict(df=1.45e-7, cf=4182, denf=998, nf=1.33),
                         "bounds": dict(thf=[-5000, 5000])},
            },
            params=[
                Param("d", 5e-7, "D amostra", "m²/s", True, 1e-8, 1e-5),
                Param("ths", 4.0, "θ amostra", "rad", True, -50, 50),
                Param("thf", 0.0, "θ fluido", "rad", True, -200, 200),
                Param("m", 22.1, "m", ""),
                Param("w0", 57e-6, "w0e", "m"),
                Param("v", 2.07, "V", ""),
                Param("l", 1.85e-3, "L espessura", "m"),
                Param("df", 2.2e-5, "D fluido", "m²/s"),
                Param("cs", 810, "c amostra", "J/(kg K)"),
                Param("cf", 1005, "c fluido", "J/(kg K)"),
                Param("dens", 2890, "ρ amostra", "kg/m³"),
                Param("denf", 1.18, "ρ fluido", "kg/m³"),
                Param("ns", 1.63, "n amostra", ""),
                Param("nf", 1.0, "n fluido", ""),
                Param("E", 72, "E (Young)", ""),
                Param("coefexp", 7.5, "α_T expansão", ""),
                Param("dndts", 5.3, "dn/dT amostra", ""),
                Param("nu", 0.23, "ν Poisson", ""),
                Param("qpen", 0.9e-3, "q⊥", ""),
                Param("qpar", 0.09e-3, "q∥", ""),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        u2 = (a * p["w0"]) ** 2
        gauss = np.exp(-u2 / 8)
        K = np.zeros(np.broadcast_shapes(t.shape, a.shape))
        # --- amostra (argfasetl) ---
        if p["ths"] != 0:
            l, ns, nf, nu, ce = p["l"], p["ns"], p["nf"], p["nu"], p["coefexp"]
            x = l * a
            hl = _hh(x) / x
            qci2 = ns ** 3 * p["E"] * ce / (4 * (1 - nu))
            B = (p["dndts"] + 4 * (ns - nf) * (1 + nu) * ce * hl
                 + qci2 * (p["qpar"] + 3 * p["qpen"]
                           - 4 * (p["qpar"] * nu + p["qpen"] * (2 + nu)) * hl))
            K = K + p["ths"] * l * gauss * B * (-np.expm1(-p["d"] * t * a * a)) / a
        # --- fluido (argfasear) ---
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
            # t0 → ±∞: os termos com t0 tendem a zero pelos dois lados (limite contínuo)
            return gauss / ((k2 - kf2) * a * a) * (
                k * erf(a * np.sqrt(df * t)) - kf * erf(a * np.sqrt(d * t)))
        t0 = (k2 - kf2) * d * df / den
        temp1 = gauss / ((k2 - kf2) * a * a)
        bb1 = k * erf(a * np.sqrt(df * t)) - kf * erf(a * np.sqrt(d * t))

        # Cada termo é  A·exp(−t0 α² t) + R.  No ramo erf do original,
        #   −c s e^{−t0α²t} erf(x) = −c s e^{−t0α²t} + c s e^{−dd α² t} erfcx(x),
        # com x² = (dd − t0) α² t. Quando os dois termos usam esse ramo,
        # k·s(df) = kf·s(d) exatamente e as partes e^{−t0α²t} se cancelam — o que
        # importa para t0 < 0 (ex.: amostra em água), onde e^{−t0α²t} estoura.
        # Ramo t0 > dd: continuação analítica do ramo erf (erf(iy)/i = erfi(y)).
        # O tllam.c original tinha este termo com o sinal trocado; conferido contra
        # a inversão numérica da solução no espaço de Laplace (tests/validate.py).
        def branch(dd, coef):
            if t0 > dd:
                return 0.0, -(2 * coef / SQRT_PI) * np.sqrt(dd / (t0 - dd)) * \
                    np.exp(-a2t * dd) * dawsn(a * np.sqrt((t0 - dd) * t))
            s = coef * np.sqrt(dd / (dd - t0))
            return -s, s * np.exp(-dd * a2t) * erfcx(a * np.sqrt((dd - t0) * t))

        A3, R3 = branch(d, kf)
        A2, R2 = branch(df, k)
        bb = R2 - R3
        if not (t0 < d and t0 < df):        # só um ramo erf: aqui t0 > 0
            bb = bb + (A2 - A3) * np.exp(-t0 * a2t)
        return temp1 * (bb1 + bb)


# ---------------------------------------------------------------------------
# TM-HAM on/off + TL no ar — equivalente a tmhamtl-onoff.c
# ---------------------------------------------------------------------------


class TMHAMTL(Model):
    linear = ("th",)     # a lente no ar também é proporcional a θ (pref ∝ th)

    def __init__(self):
        super().__init__(
            key="tmhamtl",
            title="Espelho térmico HAM on/off + TL no ar",
            description="Espelho térmico (alta absorção) com excitação desligada em "
                        "t = t_off, mais lente térmica no ar (ativa se dn/dT ar ≠ 0). "
                        "Equivalente a tmhamtl-onoff.c.",
            params=[
                Param("d", 3.46e-6, "D amostra", "m²/s", True, 1e-8, 1e-4),
                Param("th", -0.57, "θ TM", "", True, -5, 5),
                Param("m", 6.01, "m", ""),
                Param("w0", 646e-6, "w0e", "m"),
                Param("v", 52.6, "V", ""),
                Param("beta", 1.73e-5, "β expansão", "1/K"),
                Param("nu", 0.29, "ν Poisson", ""),
                Param("dndt", 0.0, "dn/dT ar", "1/K"),
                Param("df", 2.2e-5, "D ar", "m²/s"),
                Param("of", 0.899, "t_off", "s"),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
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
        """argfasear/argtau: integral em τ feita com τ = t sin²φ e tanh-sinh.

        Vetorizada em (tempo × α × φ), em blocos de `block` tempos para limitar
        a memória (~block·Nα·Nφ números por bloco)."""
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


MODELS = {m.key: m for m in (TLLAM(), TMHAMTL())}
