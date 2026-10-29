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
    # False para modelos em forma fechada (sem K·M): o ajuste em dois níveis não se aplica
    phase_based = True
    # definição física de cada θ ajustável: {param: {"def", "unit", "note"}} (mostrada na interface)
    theta_info = {}

    # parâmetros lineares cuja fase não tem a forma K·M (ex.: lente de população): a subclasse
    # devolve Φ_j (Nt, Ng) para θ_j = 1 em direct_phase(nome, t, g, p)
    direct = ()

    # implementado pelas subclasses
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
            # parâmetros fora do domínio do modelo: devolve NaN (o ajuste descarta)
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
        self.model = model
        self.m, self.w0, self.v = float(p["m"]), float(p["w0"]), float(p["v"])
        u, _ = _alpha_grid(self.m, self.w0, acc)
        tl = self.t[self.live][:, None]
        self.tl = tl
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
        direct = self.model.direct
        for n in direct:                      # termos diretos fixos (não variados)
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

    def unit_phase(self):
        """Maior |Φ_j| (em rad) por unidade de cada θ_j, na malha atual."""
        if not self.ok:
            return [1.0] * len(self.lin)
        if not hasattr(self, "base"):
            self._phase([0.0] * len(self.lin))
        return [float(np.abs(P).max()) or 1.0 for P in self.grid["Pj"]]

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
    theta_info = {
        "ths": {"def": "θ amostra = ϑ·10⁻⁶,  ϑ = P_e A_e φ / (k_s λ_p)", "unit": "10⁶ K/m",
                "note": "A fase da amostra é ϑ·l·χ(α,l), com χ (dn/dT, α_T, ν, E, q∥, q⊥) em 10⁻⁶/K. "
                        "Não é o θ em rad do modelo de Shen: a espessura l e o coeficiente χ ficam fora dele."},
        "thf": {"def": "θ fluido = 2 P_e A_e φ (∂n_f/∂T) / λ_p", "unit": "W/(K·m²)",
                "note": "θ fluido / θ amostra = 2·10⁶ k_s (∂n_f/∂T): com k_s conhecido, dá o ∂n/∂T do fluido."},
    }
    linear = ("ths", "thf")

    def __init__(self):
        super().__init__(
            key="tllam",
            title="TL+TL-Fluid (TLLAM)",
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
                Param("ths", 4.0, "θ amostra", "10⁶ K/m", True, -50, 50),
                Param("thf", 0.0, "θ fluido", "W/(K·m²)", True, -200, 200),
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
# Lente térmica de Shen (forma fechada, modo descasado)
#   Shen, Lowe & Snook, Chem. Phys. 165, 385 (1992); Baesso et al., JAP 75, 3732 (1994)
#   I(t) = I(0) [1 − (θ/2) atan(2mV / ([(1+2m)² + V²] t_c/(2t) + 1 + 2m + V²))]²
# ---------------------------------------------------------------------------


class TLShen(Model):
    theta_info = {
        "th": {"def": "θ = −P_e A_e L φ (ds/dT) / (k λ_p)", "unit": "rad",
               "note": "θ da lente térmica clássica: inclui a espessura L. Não confundir com o θ_TM do "
                       "espelho térmico (1/m) nem com o θ amostra do TLLAM. ds/dT > 0: θ < 0 e o sinal sobe."},
    }
    linear = ()
    phase_based = False     # sem fase integral: o ajuste usa PSO + LM direto
    proportional = ("th",)  # na série, θ ∝ P (inclinação θ/P)

    def __init__(self):
        super().__init__(
            key="tlshen",
            title="TL de Shen (clássico)",
            description="Expressão fechada do modelo aberrante de Shen para lente "
                        "térmica no modo descasado (amostra fina, baixa absorção, θ "
                        "pequeno). t_c = w0e²/(4D). Baesso, Shen & Snook, J. Appl. Phys. "
                        "75, 3732 (1994), eq. 10.",
            params=[
                Param("d", 1.4e-7, "D amostra", "m²/s", True, 1e-9, 1e-4),
                Param("th", -0.1, "θ", "rad", True, -3, 3),
                Param("m", 30.0, "m", ""),
                Param("w0", 50e-6, "w0e", "m"),
                Param("v", 2.0, "V", ""),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
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
# Lente térmica com fotobranqueamento — Herculano et al., J. Phys. Chem. B 117, 1932 (2013)
#   β(t) = β [(1−ε̄) e^{−K_T t} + ε̄]  (absorção que decai de β a ε̄·β com taxa K_T)
#   Φ(g,t) = θ/t_c ∫₀ᵗ [ε̄ + (1−ε̄) e^{−K_T (t−ξ)}] (1+2ξ/t_c)⁻¹ [1 − e^{−2mg/(1+2ξ/t_c)}] dξ
# Em espaço de Hankel o mesmo resultado é  Φ = ∫ K(α,t)[J0(α w0 √(m g)) − 1] dα  com
#   K = −θ e^{−α²w²/8} F(α,t)/α,
#   F = ε̄ (1 − e^{−bt}) + (1−ε̄) b (e^{−bt} − e^{−K_T t})/(K_T − b),   b = D α²
# (convolução de β(t) com a resposta térmica; K_T → 0 ou ε̄ = 1 dão a TL clássica de Shen).
# ---------------------------------------------------------------------------


def _phi1(z):
    """(1 − e^{−z})/z, estável em z → 0."""
    z = np.abs(z)
    small = z < 1e-8
    zs = np.where(small, 1.0, z)
    return np.where(small, 1.0 - 0.5 * z, -np.expm1(-zs) / zs)


class TLFotobranq(Model):
    theta_info = {
        "th": {"def": "θ = −P_e β L φ (dn/dT) / (k λ_p),  β = absorção em t = 0 (antes de branquear)",
               "unit": "rad",
               "note": "Mesmo θ da TL de Shen, mas com o β inicial: o sinal em t longo corresponde a ε̄·θ. "
                       "dn/dT < 0 (água): θ > 0 e o sinal cai. Na série, θ ∝ P_e. "
                       "K_T [1/s] = K_m + K_pb P_e (K_m: difusão de massa; K_pb: fotobranqueamento) "
                       "e ε̄ = ε_produtos/ε_reagentes — a reta de K_T × P_e aparece na série."},
    }
    linear = ("th",)

    def __init__(self):
        super().__init__(
            key="tlbranq",
            title="TL com fotobranqueamento (Herculano 2013)",
            description="Lente térmica de líquido com absorção que decai por reação fotoquímica "
                        "(pseudo-primeira ordem, taxa K_T) até ε̄·β. Ajusta θ, K_T e ε̄; fixe D com uma "
                        "medida do solvente puro (K_T = 0). Herculano et al., J. Phys. Chem. B 117, "
                        "1932 (2013), eqs. 3–6.",
            params=[
                Param("d", 1.41e-7, "D amostra", "m²/s", False, 1e-9, 1e-4),
                Param("th", 0.05, "θ", "rad", True, -5, 5),
                Param("kt", 2.0, "K_T", "1/s", True, 1e-4, 1e3),
                Param("eb", 0.3, "ε̄", "", True, 0, 5),
                Param("m", 35.1, "m", ""),
                Param("w0", 70e-6, "w0e", "m"),
                Param("v", 7.5, "V", ""),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
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
# Lente térmica + lente de população — Silva et al., Appl. Phys. Lett. 102, 141910 (2013)
#   Φ = Φ_th + Φ_PL
#   Φ_th = θ_th ∫₀ᵗ dξ/(t_c + 2ξ) [1 − e^{−2mg/(1+2ξ/t_c)}]            (= TL de Shen)
#   Φ_PL = θ_PL ∫₀ˡ  u/(1/τ + h u) [1 − e^{−t(1/τ + h u)}] dz,   u = e^{−2mg} e^{−A_e z}
# (população do estado excitado por equação de taxa de dois níveis, saturável; h é a taxa de
# excitação no centro, h = I0 (σ_abs + σ_se)/E). Com u = e^{−2mg − A_e z} a integral em z fica
#   Φ_PL = θ_PL /(A_e h) ∫_{s_L}^{s_0} (1 − e^{−t s})/s ds,  s = 1/τ + h u,
#   ∫ = ln(s0/sL) − [E1(t sL) − E1(t s0)]   (série em t s para t s0 ≤ 1; quadratura em z se A_e l < 0,02).
# Como no artigo, a fase é Φ_PL(g) − Φ_PL(0) (uma fase constante não altera a intensidade).
# ---------------------------------------------------------------------------

from scipy.special import exp1 as _exp1

_GL8 = np.polynomial.legendre.leggauss(8)


def _pl_unit_phase(t, u0, ae, l, h, tau):
    """Φ_PL para θ_PL = 1; t (Nt,1), u0 = e^{−2mg} (1,Ng) ou (Ng,)."""
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
    # série: Σ (−1)^{j+1} t^j (s0^j − sL^j)/(j·j!)
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
               "note": "Mesmo θ do Shen. ds/dT < 0 (ZBLAN) com aquecimento: θ_th > 0; "
                       "resfriamento a laser inverte o sinal (φ < 0)."},
        "tpl": {"def": "θ_PL = I0 σ_abs N_T Δα_P (2π f_L)² / (n λ_p E)",
                "unit": "1/(m·s)",
                "note": "Força da lente de população. Δα_P: diferença de polarizabilidade eletrônica "
                        "excitado–fundamental; f_L = (n²+2)/3; E = h c/λ_ex; I0 = 2P/(π w_e²)."},
    }
    linear = ("th", "tpl")
    direct = ("tpl",)

    def __init__(self):
        super().__init__(
            key="tlpop",
            title="TL + lente de população (Silva 2013)",
            description="Lente térmica mais lente de população do estado excitado (vida τ, taxa de "
                        "excitação h = I0(σ_abs+σ_se)/E no centro, absorção A_e, espessura l). "
                        "Silva et al., Appl. Phys. Lett. 102, 141910 (2013), eqs. 1–3. "
                        "Ajusta θ_th, θ_PL, τ e D; m e V são da montagem, h, A_e e l são conhecidos.",
            params=[
                Param("d", 2.5e-7, "D amostra", "m²/s", True, 1e-9, 1e-4),
                Param("th", 0.04, "θ_th", "rad", True, -3, 3),
                Param("tpl", 2900.0, "θ_PL", "1/(m·s)", True, -1e7, 1e7),
                Param("tau", 1.7e-3, "τ (vida)", "s", True, 1e-5, 1e-1),
                Param("h", 40.0, "h (taxa de excitação)", "1/s", False, 1e-3, 1e5),
                Param("ae", 240.0, "A_e absorção", "1/m", False, 0.0, 1e5),
                Param("l", 3.79e-3, "l espessura", "m"),
                Param("m", 28.7, "m", ""),
                Param("w0", 144e-6, "w0e", "m"),
                Param("v", 15.0, "V", ""),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
            ],
        )

    def kernel(self, t, a, p, acc):
        # parte térmica (Shen): K = −θ_th e^{−α²w²/8} (1 − e^{−Dα²t})/α
        return -p["th"] * np.exp(-(a * p["w0"]) ** 2 / 8) * (-np.expm1(-p["d"] * a * a * t)) / a

    def direct_phase(self, name, t, g, p):
        g = np.asarray(g, float)
        u0 = np.exp(-2 * p["m"] * g)
        both = np.concatenate([[1.0], u0])          # u = 1 em g = 0 (referência do centro)
        ph = _pl_unit_phase(t, both, p["ae"], p["l"], p["h"], p["tau"])
        return ph[:, 1:] - ph[:, :1]


# ---------------------------------------------------------------------------
# Espelho térmico: LAM, BAM (lei de Beer) e LAM com espessura finita
#   Sato et al., J. Appl. Phys. 104, 053520 (2008)  — LAM, BAM, HAM
#   Astrath et al., J. Opt. Soc. Am. B 28, 1735 (2011) — espessura finita
#
# Com θ_TM = −P_e A_e α_T (1+ν) φ/(λ_p k) [1/m], a fase é
#   Φ(g,t) = 4 θ_TM ∫ e^{−α²w²/8} F(α,t) [J0(α w √(mg)) − 1] dα,
#   F_LAM  = (1 − e^{−Dα²t})/α²                       (Sato eq. 24/32)
#   F_fin  = h(lα) (1 − e^{−Dα²t})/α²,  h(x) = (cosh x − 1)/(x + sinh x)
#                                                      (Astrath 2011, eq. 10)
#   F_BAM  = α D Ĝ(α,t),  Ĝ = ∫_0^∞ T̂(α,z,t) e^{−αz} dz  (fonte e^{−A_e z})
# Ĝ sai de dĜ/dt = [A_e erfc(αs) − α e^{−α²s²} erfcx(A_e s)]/(A_e² − α²),
# s = √(Dt), integrada em forma fechada (ou por série quando s·max(α, A_e) ≤ 0,75).
# Limites: A_e → 0 recupera o LAM; A_e → ∞ recupera o HAM (θ_HAM = θ_TM/A_e).
# Desligamento em t_off por superposição: F(t) − F(t − t_off).
# ---------------------------------------------------------------------------

from math import factorial as _factorial
from scipy.special import gamma as _gamma_fn

_XERFC_C = np.array([(-1) ** n / (_factorial(n) * (2 * n + 1) * (2 * n + 3)) for n in range(14)])


def _xerfc_int(X):
    """∫_0^X x erfc(x) dx (série de Taylor para X pequeno, sem cancelamento)."""
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
    """Ĝ = t/(A_e+α) Σ_n p_n ŝ^n/(n/2+1), em variáveis escaladas por σ = max(α, A_e).

    Vem de erfcx(z) = Σ (−z)^n/Γ(n/2+1): o fator (A_e² − α²) cancela termo a termo,
    então não há divisão por números pequenos perto de α = A_e. Os coeficientes p_n
    dependem só de α: são calculados uma vez por α distinto."""
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
    for k in range(0, n, 2):                 # E só tem termos pares
        P[k:] += E[k] * A[:n - k]
    P /= (np.arange(n) / 2 + 1)[:, None]
    P = P[:, inv]
    sh = sig[inv] * np.sqrt(d * t)
    acc = np.zeros(a.size)
    for k in range(n - 1, -1, -1):
        acc = acc * sh + P[k]
    return t / (ae + a) * acc


def bam_g(a, t, ae, d, rel=1e-4, zmax=0.75):
    """Ĝ(α,t) do modelo BAM (fonte e^{−A_e z}, superfície livre e adiabática).

    Conferido contra a quadratura direta de dĜ/dt (erro relativo < 3e-7 em
    10⁻¹ < A_e < 10⁸ m⁻¹, 10⁻⁶ < t < 10 s, inclusive α ≈ A_e)."""
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
        if np.any(near):   # singularidade removível em α = A_e: interpolação linear
            tn = tcl[near]
            lo = _bam_closed(np.full(tn.size, ae * (1 - rel)), tn, ae, d)
            hi = _bam_closed(np.full(tn.size, ae * (1 + rel)), tn, ae, d)
            g[near] = 0.5 * (lo + hi) + (ac[near] - ae) / (2 * rel * ae) * (hi - lo)
        out[cl] = g
    return out.reshape(shape)


class _TMBase(Model):
    """Espelho térmico (fase 4 θ_TM ∫ e^{−α²w²/8} F(α,t) [J0 − 1] dα) com on/off."""
    linear = ("th",)
    theta_info = {
        "th": {"def": "θ_TM = −P_e A_e α_T (1+ν) φ / (λ_p k)", "unit": "1/m",
               "note": "Mesmo θ_TM nos modelos LAM, BAM e amostra finita (é o θ_TM dos artigos, em m⁻¹). "
                       "O modelo HAM usa θ_TM/A_e, adimensional. Amostra que dilata (α_T > 0): θ_TM < 0."},
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
    Param("d", 6e-7, "D amostra", "m²/s", True, 1e-8, 1e-4),
    Param("th", -300.0, "θ TM", "1/m", True, th_lo, th_hi),
    Param("m", 20.0, "m", ""),
    Param("w0", 44e-6, "w0e", "m"),
    Param("v", 5.2, "V", ""),
]


class TMLAM(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmlam",
            title="TM LAM (baixa absorção)",
            description="Espelho térmico para amostra semi-infinita com fonte uniforme "
                        "em z (A_e·espessura ≪ 1). LAM = Low Absorption Model. Sato et al., J. Appl. Phys. 104, "
                        "053520 (2008), eq. 32. θ_TM = −P_e A_e α_T(1+ν)φ/(λ_p k), em 1/m.",
            params=_TM_COMMON(-1e6, 1e6) + [
                Param("of", 0.0, "t_off (0 = sem desligar)", "s"),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        return -np.expm1(-p["d"] * a * a * t) / (a * a)


class TMFinito(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmfinito",
            title="TM LAM, amostra finita",
            description="Espelho térmico de baixa absorção em amostra de espessura l, "
                        "com as duas faces livres. Astrath et al., J. Opt. Soc. Am. B 28, "
                        "1735 (2011), eq. 10. Tende ao LAM semi-infinito para l grande.",
            params=_TM_COMMON(-1e6, 1e6) + [
                Param("l", 1e-3, "l espessura", "m"),
                Param("of", 0.0, "t_off (0 = sem desligar)", "s"),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        return _hh(p["l"] * a) * -np.expm1(-p["d"] * a * a * t) / (a * a)


class TMBAM(_TMBase):
    def __init__(self):
        super().__init__(
            key="tmbam",
            title="TM BAM (lei de Beer)",
            description="Espelho térmico com absorção pela lei de Beer, válido para "
                        "qualquer A_e (amostra semi-infinita); BAM = Beer Absorption Model. Sato et al., J. Appl. Phys. "
                        "104, 053520 (2008), eq. 33–34. Tende ao LAM (A_e → 0) e ao HAM "
                        "(A_e → ∞, com θ_HAM = θ_TM/A_e).",
            params=_TM_COMMON(-1e9, 1e9) + [
                Param("ae", 1000.0, "A_e absorção", "1/m", False, 1.0, 1e8),
                Param("of", 0.0, "t_off (0 = sem desligar)", "s"),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
            ],
        )

    def _F(self, t, a, p):
        d = p["d"]
        return a * d * bam_g(a, t, float(p["ae"]), d)


# ---------------------------------------------------------------------------
# TM+TL-Fluid, HAM on/off (aproximação, sem retorno do fluido) — equivalente a tmhamtl-onoff.c
# ---------------------------------------------------------------------------


class TMHAMTL(Model):
    theta_info = {
        "th": {"def": "θ_HAM = θ_TM / A_e = −P_e α_T (1+ν) φ / (λ_p k)", "unit": "adimensional",
               "note": "Não contém A_e (amostra opaca). Para comparar com LAM, BAM e amostra finita, "
                       "que dão θ_TM em 1/m: θ_TM = θ_HAM·A_e. Amostra que dilata (α_T > 0): θ < 0. "
                       "Com P_e incidente e refletividade R, o θ_TM de Flizikowski et al. (2020) é (1−R)·θ_HAM."},
    }
    linear = ("th",)     # a lente no ar também é proporcional a θ (pref ∝ th)

    def __init__(self):
        super().__init__(
            key="tmhamtl",
            title="TM+TL-Fluid, HAM on/off (aproximação: sem retorno do fluido)",
            description="Espelho térmico HAM (High Absorption Model, alta absorção) com excitação desligada em "
                        "t = t_off, mais a lente térmica do fluido (ar), ativa se dn/dT ≠ 0. "
                        "Equivalente a tmhamtl-onoff.c: aproximação em que a superfície da amostra é "
                        "livre; o fluido recebe calor e forma a TL, mas não esfria a amostra de volta. "
                        "Para o acoplamento completo (Flizikowski et al., Opt. Express 28, 7116, 2020) "
                        "use TM+TL-Fluid, HAM (acoplamento completo): nos três metais do artigo a "
                        "diferença é de 1e-4 a 1e-3 do sinal com o laser ligado e chega a ~1% logo "
                        "após desligar (inox em ar).",
            params=[
                Param("d", 3.46e-6, "D amostra", "m²/s", True, 1e-8, 1e-4),
                Param("th", -0.57, "θ HAM", "adim.", True, -5, 5),
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


# ---------------------------------------------------------------------------
# Espelho térmico (amostra opaca) com lente térmica no fluido acoplado — com acoplamento
# de calor amostra–fluido completo.
#   Lukasievicz et al., Appl. Spectrosc. 66, 1461 (2012), eqs. 4–9, 18–22
#
# Fonte superficial; amostra em z>0 e fluido em z<0. Em Hankel–Laplace (q_i² = α² + s/D_i)
#   T(α, z=0, s) ∝ 1 / [s (k_s q_s + k_f q_f)],
#   Ĝ_s = ∫_0^∞ T e^{−αz} dz ∝ 1/[s (k_s q_s + k_f q_f)(q_s + α)]   (espelho)
#   Ĝ_f = ∫_{−∞}^0 T dz     ∝ 1/[s q_f (k_s q_s + k_f q_f)]         (lente no fluido)
# Com ŝ = s/(D_s α²), τ = D_s α² t, κ = k_f/k_s e r = D_s/D_f, ambos dependem só de τ:
#   g_s(τ) = L⁻¹{1/[ŝ (√(1+ŝ) + κ√(1+rŝ)) (√(1+ŝ) + 1)]}
#   g_f(τ) = L⁻¹{1/[ŝ √(1+rŝ) (√(1+ŝ) + κ√(1+rŝ))]}
# tabelados em log τ por inversão numérica de Laplace (Talbot) e interpolados.
# Fase  Φ = ∫ K [J0(α w0 √(mg)) − 1] dα,  K = e^{−α²w²/8} (4 h_TM g_s + (h_f/k_s) g_f)/α.
# Sem fluido (κ = 0) K_TM é o do HAM (θ_HAM = h_TM).
# ---------------------------------------------------------------------------

from scipy.interpolate import CubicSpline as _CubicSpline

_TAU_LO, _TAU_HI, _TAU_N = -12.0, 9.5, 2600
_TM_TABLES: dict = {}


def _talbot(F, tau, M=22):
    """Inversa de Laplace (Abate–Valkó, Talbot fixo) de F(ŝ) em vários τ; ~1e-12."""
    r = 2 * M / (5 * tau)
    th = np.arange(1, M) * np.pi / M
    cot = 1 / np.tan(th)
    s = r[None, :] * th[:, None] * (cot[:, None] + 1j)
    sig = th + (th * cot - 1) * cot
    f0 = 0.5 * np.real(F(r.astype(complex))) * np.exp(r * tau)
    terms = np.real(np.exp(tau[None, :] * s) * F(s) * (1 + 1j * sig[:, None]))
    return (r / M) * (f0 + terms.sum(0))


def _tm_table(kind, kappa, rr):
    """Spline de g_s ('s') ou g_f ('f') em log10 τ, para (κ, r) dados."""
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
    """g_s ou g_f em τ qualquer (assíntotas fora da tabela)."""
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
        "th": {"def": "h_TM = θ_HAM = −P_e α_T (1+ν) φ / (λ_p k_s)", "unit": "adimensional",
               "note": "Mesmo θ do TM-HAM (amostra opaca): não contém A_e. Amostra que dilata: h_TM < 0. "
                       "Com P_e incidente e refletividade R, o θ_TM de Flizikowski et al. (2020) é (1−R)·h_TM."},
        "thf": {"def": "h_f = 2 P_e φ (∂n/∂T)_f / λ_p", "unit": "W/(K·m)",
                "note": "Lente térmica do fluido (dois passes: antes e depois da reflexão). "
                        "Com k_s conhecido, h_f/P_e dá (∂n/∂T)_f: h_f = 2 P_e φ (∂n/∂T)_f/λ_p. "
                        "P_e é a potência absorvida; com P_e incidente use P_e(1−R)."},
    }
    linear = ("th", "thf")

    def __init__(self):
        super().__init__(
            key="tmfluido",
            title="TM+TL-Fluid, HAM (acoplamento completo, Flizikowski 2020)",
            description="Espelho térmico de amostra opaca e lente térmica do fluido não absorvente "
                        "em contato, com troca de calor amostra–fluido (sem a restrição k_f√D ≪ k√D_f). "
                        "Lukasievicz et al., Appl. Spectrosc. 66, 1461 (2012), eqs. 20–22. "
                        "k_s = D_s ρ_s c_s e k_f = D_f ρ_f c_f: fixe ρ e c e ajuste as difusividades. "
                        "Para ar, o fluido quase não muda a amostra; em água, sim. "
                        "É o modelo completo de Flizikowski et al., Opt. Express 28, 7116 (2020), eqs. 3–10; "
                        "o TM+TL-Fluid, HAM on/off (aproximação) ignora o retorno do fluido sobre a amostra.",
            presets={
                "Ar": {"values": dict(df=2.19e-5, cf=1005, denf=1.18),
                       "bounds": dict(thf=[-1e3, 1e3], df=[1e-6, 1e-3])},
                "Água": {"values": dict(df=1.45e-7, cf=4182, denf=998),
                         "bounds": dict(thf=[-1e5, 1e5], df=[1e-8, 1e-5])},
            },
            params=[
                Param("d", 6.94e-6, "D amostra", "m²/s", True, 1e-8, 1e-4),
                Param("th", -0.5, "h_TM (θ HAM)", "adim.", True, -50, 50),
                Param("thf", -100.0, "h_f", "W/(K·m)", True, -1e5, 1e5),
                Param("df", 1.45e-7, "D fluido", "m²/s", False, 1e-8, 1e-3),
                Param("m", 32.8, "m", ""),
                Param("w0", 74e-6, "w0e", "m"),
                Param("v", 6.9, "V", ""),
                Param("cs", 460, "c amostra", "J/(kg K)"),
                Param("dens", 7800, "ρ amostra", "kg/m³"),
                Param("cf", 4182, "c fluido", "J/(kg K)"),
                Param("denf", 998, "ρ fluido", "kg/m³"),
                Param("of", 0.0, "t_off (0 = sem desligar)", "s"),
                Param("scale", 1.0, "I0 (escala)", "", False, 0.8, 1.2),
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



# Ordem na interface: primeiro os modelos de TL, depois os de TM (o primeiro é o padrão ao abrir).
MODELS = {m.key: m for m in (TLLAM(), TLShen(), TLFotobranq(), TLPopulacao(),
                             TMLAM(), TMFinito(), TMBAM(), TMHAMTL(), TMFluido())}
