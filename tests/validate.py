"""Validação do motor novo.

1. Compara com curvas geradas pelos programas originais (tllam.c e
   tmhamtl-onoff.c, compilados com GSL) — referências em tests/ref/.
   O original imprime 6 casas decimais e integra com tolerância relativa
   1e-4/1e-5, então diferenças de até ~1e-5 são esperadas.
2. Autoconvergência: resultado com quadraturas 2× mais finas.

Uso: python tests/validate.py
"""
import sys
from pathlib import Path

import numpy as np

# símbolos como − e χ² não existem na codificação antiga do Windows (cp1252);
# sem isto, imprimir o resultado aborta o teste quando a saída é redirecionada
for _s in (sys.stdout, sys.stderr):
    if hasattr(_s, "reconfigure"):
        _s.reconfigure(errors="replace")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tlfit.models import MODELS, Accuracy  # noqa: E402

REF = ROOT / "tests" / "ref"
CASES = [
    ("tllam", "tllam_d5.317e-7_ths4.33_thf0.dat", dict(d=5.317228739232711e-07, ths=4.330896764158691, thf=0)),
    ("tllam", "tllam_d5e-7_ths4_thf30.dat", dict(d=5e-7, ths=4, thf=30)),
    ("tmhamtl", "tmhamtl_default.dat", dict()),
]
ok = True
for key, ref, upd in CASES:
    m = MODELS[key]
    p = m.defaults()
    p.update(upd)
    r = np.loadtxt(REF / ref, encoding="utf-8")
    t = np.ascontiguousarray(r[:, 0])
    I = m.intensity(t, p)
    I2 = m.intensity(t, p, Accuracy().scaled(2))
    d_orig, d_self = np.abs(I - r[:, 1]).max(), np.abs(I - I2).max()
    print(f"{key:8s} {ref:34s} |novo−original| {d_orig:.1e}   autoconvergência {d_self:.1e}")
    ok &= d_orig < 2e-4 and d_self < 1e-7

# 4. Termo do fluido (TLLAM) contra a inversão numérica (Talbot) da solução exata
#    no espaço de Laplace: amostra e fluido semi-infinitos, fonte uniforme na amostra,
#    ∫T_f dz ∝ 1/(s q_s q_f (k q_s + k_f q_f)),  q² = α² + s/D.
#    Cobre os dois ramos da fórmula fechada: Dawson quando (e_s − e_f)(D − D_f) > 0.
def talbot(F, t, M=32):
    r = 2 * M / (5 * t)
    th = np.arange(1, M) * np.pi / M
    cot = 1 / np.tan(th)
    s = r * th * (cot + 1j)
    sig = th + (th * cot - 1) * cot
    return (r / M) * (0.5 * np.real(F(r)) * np.exp(r * t)
                      + np.sum(np.real(np.exp(t * s) * F(s) * (1 + 1j * sig))))


rng = np.random.default_rng(1)
worst = {"ramo erf": 0.0, "ramo Dawson": 0.0}
for _ in range(200):
    p = dict(cs=rng.uniform(300, 3000), dens=rng.uniform(800, 8000), cf=rng.uniform(800, 4500),
             denf=rng.choice([1.2, 800, 1000, 1300]), d=10 ** rng.uniform(-8, -5),
             df=10 ** rng.uniform(-8, -4.5))
    k, kf = p["d"] * p["cs"] * p["dens"], p["df"] * p["cf"] * p["denf"]
    t, a = 10 ** rng.uniform(-4, 0), rng.uniform(0.05, 6) / 57e-6
    qs = lambda s: np.sqrt(a * a + s / p["d"])       # noqa: E731
    qf = lambda s: np.sqrt(a * a + s / p["df"])      # noqa: E731
    ex = a * talbot(lambda s: 1 / (s * qs(s) * qf(s) * (k * qs(s) + kf * qf(s))), t)
    v = MODELS["tllam"]._fluid(np.array([[t]]), np.array([[a]]), p, 1.0)[0, 0]
    es, ef = k / np.sqrt(p["d"]), kf / np.sqrt(p["df"])        # efusividades
    key = "ramo Dawson" if (es - ef) * (p["d"] - p["df"]) > 0 else "ramo erf"
    worst[key] = max(worst[key], abs(v - ex) / abs(ex))
for key, e in worst.items():
    print(f"fluido   Laplace, {key}: erro relativo máx {e:.1e}")
    ok &= e < 1e-6

# 5. Ajuste completo (modo automático, dois níveis) da amostra em água:
#    deve achar sempre o mesmo mínimo, sem depender de sorte.
from tlfit.optim import FitProblem, make_free, run_fit  # noqa: E402

dat = np.loadtxt(ROOT / "tlfit" / "examples" / "tllam_amostra_agua.dat", encoding="utf-8")
m = MODELS["tllam"]
p = m.defaults()
p.update(df=1.45e-7, cf=4182, denf=998, nf=1.33, d=5e-7, ths=4, thf=0)
prob = FitProblem(m, dat[:, 0], dat[:, 1], p, [make_free("d", 1e-8, 1e-5),
                  make_free("ths", -50, 50), make_free("thf", -5000, 5000)])
r = run_fit(prob, method="auto", log=lambda s: None)
print(f"ajuste   água, modo automático: χ² = {r['chi2']:.5e}, D = {r['params']['d']:.5e}, "
      f"θs = {r['params']['ths']:.4f}, θf = {r['params']['thf']:.1f} ({r['seconds']:.1f} s)")
ok &= abs(r["chi2"] / 9.0531e-6 - 1) < 1e-3 and abs(r["params"]["d"] / 4.9774e-7 - 1) < 1e-3

# 6. Série sintética (vidro em ar, θ_amostra/P = 0,020 rad/mW, D = 5,4e-7):
#    a inclinação pela origem deve recuperar θ/P e D não deve ter tendência com P.
import glob  # noqa: E402
import os  # noqa: E402

from tlfit.batch import fit_series, series_stats  # noqa: E402
from tlfit.dataio import guess_power, parse_table  # noqa: E402

items = []
for fn in sorted(glob.glob(str(ROOT / "tlfit" / "examples" / "serie_vidro_ar" / "*.dat"))):
    txt = open(fn, encoding="utf-8").read()
    d, _ = parse_table(txt)
    items.append(dict(name=os.path.basename(fn), t=d[:, 0], y=d[:, 1],
                      x=guess_power(os.path.basename(fn), txt)[0]))
p = m.defaults()
p.update(d=5e-7, ths=4, thf=0)
free = [dict(name="d", lo=1e-8, hi=1e-5), dict(name="ths", lo=-50, hi=50),
        dict(name="thf", lo=-200, hi=200)]
st = series_stats(m, fit_series(m, items, p, free, log=lambda s: None))
a, ea = st["ths"]["slope0"], st["ths"]["slope0_err"]
print(f"série    θs/P = {a:.5f} ± {ea:.5f} rad/mW (verdadeiro 0,020); "
      f"D = {st['d']['mean']:.4e} ± {st['d']['std']:.1e} (verdadeiro 5,4e-7)")
ok &= abs(a - 0.020) < 3 * ea and abs(st["d"]["mean"] - 5.4e-7) < 3 * st["d"]["std"]

# 7. Espelho térmico BAM: Ĝ(α,t) em forma fechada/série contra a quadratura direta
#    de dĜ/dt = [A_e erfc(αs) − α e^{−α²s²} erfcx(A_e s)]/(A_e² − α²), s = √(Dτ).
from scipy.integrate import quad  # noqa: E402
from scipy.special import erfc as _erfc, erfcx as _erfcx, i0e, i1e, j0  # noqa: E402

from tlfit.models import bam_g  # noqa: E402

xg, wg = np.polynomial.legendre.leggauss(40)
edges = np.concatenate([[0], np.logspace(-8, 0, 60)])


def g_ref(a, t, ae, d):
    tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):            # τ = t x², x ∈ [0, 1]
        x = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
        s = np.sqrt(d * t) * x
        f = (ae * _erfc(a * s) - a * np.exp(-a * a * s * s) * _erfcx(ae * s)) / (ae * ae - a * a)
        tot += np.sum(0.5 * (hi - lo) * wg * f * 2 * t * x)
    return tot


rng = np.random.default_rng(2)
worst = 0.0
for _ in range(400):
    ae, t = 10 ** rng.uniform(-1, 8), 10 ** rng.uniform(-6, 1)
    a = 10 ** rng.uniform(-1, 6.5) if rng.random() < 0.7 else ae * (1 + 10 ** rng.uniform(-5, -0.5) * rng.choice([-1, 1]))
    v = bam_g(np.array([a]), np.array([t]), ae, 6e-7)[0]
    worst = max(worst, abs(v / g_ref(a, t, ae, 6e-7) - 1))
print(f"TM-BAM   Ĝ(α,t) contra quadratura em τ: erro relativo máx {worst:.1e}")
ok &= worst < 1e-6

# 8. TM-LAM: fase pelo núcleo em α contra a forma fechada em g (Sato 2008, eq. 32)
D_, w_, m_ = 6e-7, 44e-6, 10.0
tc_ = w_ * w_ / (4 * D_)


def eq32(g, t):
    x, x0 = g * m_ / (1 + 2 * t / tc_), g * m_
    A = (1 + 2 * t / tc_) ** -0.5 * ((1 + 2 * g * m_ + 2 * t / tc_) * i0e(x) + 2 * g * m_ * i1e(x))
    B = (1 + 2 * g * m_) * i0e(x0) + 2 * g * m_ * i1e(x0)
    return np.sqrt(2 * np.pi) * w_ * (A - B)


worst = 0.0
for t in (0.5 * tc_, 5 * tc_, 50 * tc_):
    for g in (0.3, 1.0, 3.0):
        K = lambda a: MODELS["tmlam"].kernel(np.array(t), a, dict(d=D_, th=1.0, w0=w_), None)  # noqa: E731
        num = quad(lambda a: K(a) * (j0(a * w_ * np.sqrt(m_ * g)) - 1), 0, 40 / w_,
                   limit=800, epsrel=1e-11, epsabs=1e-14)[0]
        ref = eq32(g, t) - eq32(0, t)
        worst = max(worst, abs(num / ref - 1))
print(f"TM-LAM   fase contra a eq. 32 de Sato et al. (2008): erro relativo máx {worst:.1e}")
ok &= worst < 1e-7

# 9. Limites: BAM → HAM (A_e grande, θ_HAM = θ_TM/A_e, com desligamento) e BAM → LAM
#    (A_e → 0); amostra finita → LAM semi-infinito (l grande).
t = np.linspace(0, 1.5, 300)
ham = MODELS["tmhamtl"].defaults()
I_ham = MODELS["tmhamtl"].intensity(t, ham)
geo = {k: ham[k] for k in ("d", "m", "w0", "v", "of", "scale")}
I_bam = MODELS["tmbam"].intensity(t, dict(geo, th=ham["th"] * 1e8, ae=1e8))
e_ham = np.abs(I_bam - I_ham).max() / np.abs(I_ham - 1).max()
t = np.linspace(0, 0.2, 300)
geo = dict(d=6e-7, th=-120.0, m=71.0, w0=44e-6, v=5.2, of=0.0, scale=1.0)
I_lam = MODELS["tmlam"].intensity(t, geo)
e_lam = np.abs(MODELS["tmbam"].intensity(t, dict(geo, ae=1e-4)) - I_lam).max() / np.abs(I_lam - 1).max()
e_fin = np.abs(MODELS["tmfinito"].intensity(t, dict(geo, l=1.0)) - I_lam).max() / np.abs(I_lam - 1).max()
print(f"limites  BAM→HAM {e_ham:.1e}   BAM→LAM {e_lam:.1e}   finito→LAM {e_fin:.1e} (relativo ao sinal)")
ok &= e_ham < 1e-7 and e_lam < 1e-5 and e_fin < 1e-9

# 10. Shen (forma fechada) contra a integral de difração com a mesma fase, θ pequeno
from tlfit.models import Model  # noqa: E402


class _TLNum(Model):
    def kernel(self, t, a, p, acc):
        return -p["th"] * np.exp(-(a * p["w0"]) ** 2 / 8) * (-np.expm1(-p["d"] * a * a * t)) / a


t = np.linspace(0, 0.5, 300)
ps = dict(d=1.4e-7, th=0.01, m=30.0, w0=50e-6, v=2.0, scale=1.0)
I_s = MODELS["tlshen"].intensity(t, ps)
e_s = np.abs(_TLNum("x", "", "").intensity(t, ps) - I_s).max() / np.abs(I_s - 1).max()
print(f"TL Shen  contra a integral de difração (θ = 0,01): diferença {e_s:.1e} do sinal (ordem θ)")
ok &= e_s < 2e-3

# 11. Ajustes dos exemplos sintéticos dos modelos novos (partindo dos valores padrão)
import re  # noqa: E402

for fn in ("tmlam_vidro_sintetico.dat", "tmbam_vidro_sintetico.dat",
           "tmfinito_vidro_sintetico.dat", "tlshen_vidro_sintetico.dat"):
    txt = (ROOT / "tlfit" / "examples" / fn).read_text(encoding="utf-8")
    key = re.search(r"# tlfit: modelo=(\w+)", txt).group(1)
    true = {k: float(v) for k, v in re.findall(r"\b(d|th) = (\S+?)[;,]", txt)}
    mm = MODELS[key]
    p = mm.defaults()
    p.update({k: float(v) for k, v in re.findall(r"(\w+)=(\S+)", txt.splitlines()[1].split("modelo=")[1])
              if k != key})
    dd, _ = parse_table(txt)
    thp = next(q for q in mm.params if q.name == "th")
    prob = FitProblem(mm, dd[:, 0], dd[:, 1], p, [make_free("d", 1e-8, 1e-4), make_free("th", thp.lo, thp.hi)])
    r = run_fit(prob, method="auto", log=lambda s: None)
    zd = (r["params"]["d"] - true["d"]) / r["sigma"]["d"]
    zt = (r["params"]["th"] - true["th"]) / r["sigma"]["th"]
    print(f"ajuste   {key:9s} D = {r['params']['d']:.4e} ({zd:+.1f}σ)  θ = {r['params']['th']:.5g} "
          f"({zt:+.1f}σ)  ({r['seconds']:.1f} s)")
    ok &= abs(zd) < 4 and abs(zt) < 4

# 12. TL com fotobranqueamento (Herculano 2013): fase em α contra a eq. 4 (quadratura em ξ),
#     limites K_T → 0 e ε̄ = 1 (Shen) e ajuste do exemplo sintético
from tlfit.models import _alpha_grid, _bessel_matrix, g_base_edges  # noqa: E402

mb = MODELS["tlbranq"]
pb = mb.defaults()
accb = Accuracy()
u_b, _ = _alpha_grid(pb["m"], pb["w0"], accb)
tb = np.array([1e-3, 5e-3, 2e-2, 0.1, 0.5, 2.0])
Kb = mb.kernel(tb[:, None], (u_b / pb["w0"])[None, :], pb, accb)
Mb, gb, _ = _bessel_matrix(pb["m"], pb["w0"], accb, pb["v"], (1,) * (len(g_base_edges(pb["v"], accb)) - 1))
Pb = Kb @ Mb
tc_b = pb["w0"] ** 2 / (4 * pb["d"])


def _eq4(g, tt):
    f = lambda x: (pb["eb"] + (1 - pb["eb"]) * np.exp(-pb["kt"] * (tt - x))) / (1 + 2 * x / tc_b) \
        * (1 - np.exp(-2 * pb["m"] * g / (1 + 2 * x / tc_b)))
    return pb["th"] / tc_b * quad(f, 0, tt, limit=200, epsabs=1e-13, epsrel=1e-12)[0]


worst = 0.0
for i, tt in enumerate(tb):
    for j in range(0, len(gb), max(1, len(gb) // 12)):
        ref = _eq4(gb[j], tt)
        worst = max(worst, abs(Pb[i, j] - ref) / max(abs(ref), 1e-3 * pb["th"]))
print(f"TL branq fase contra a eq. 4 de Herculano et al. (2013): erro relativo máx {worst:.1e}")
ok &= worst < 1e-8

sh_p = dict(d=pb["d"], th=0.01, m=pb["m"], w0=pb["w0"], v=pb["v"])
I_sh = MODELS["tlshen"].intensity(tb, sh_p)
e_k0 = np.abs(mb.intensity(tb, dict(pb, th=0.01, kt=1e-9)) - I_sh).max() / np.abs(I_sh - 1).max()
e_e1 = np.abs(mb.intensity(tb, dict(pb, th=0.01, eb=1.0)) - I_sh).max() / np.abs(I_sh - 1).max()
print(f"limites  K_T→0 {e_k0:.1e}   ε̄=1 {e_e1:.1e} (relativo ao sinal; Shen é de ordem θ)")
ok &= e_k0 < 2e-3 and e_e1 < 2e-3

txt = (ROOT / "tlfit" / "examples" / "tlbranq_eosina_sintetico.dat").read_text(encoding="utf-8")
dd, _ = parse_table(txt)
pp = mb.defaults()
pp.update(d=1.41e-7, m=35.1, w0=7e-5, v=7.5)
true_b = dict(th=0.048, kt=2.0, eb=0.27)
prob = FitProblem(mb, dd[:, 0], dd[:, 1], pp, [make_free("th", -5, 5), make_free("kt", 1e-4, 1e3), make_free("eb", 0, 5)])
r = run_fit(prob, method="auto", log=lambda s: None)
zs = {k: (r["params"][k] - true_b[k]) / r["sigma"][k] for k in true_b}
print("ajuste   tlbranq   " + "  ".join(f"{k} = {r['params'][k]:.5g} ({zs[k]:+.1f}σ)" for k in true_b)
      + f"  ({r['seconds']:.1f} s)")
ok &= all(abs(z) < 4 for z in zs.values())

# 13. TM + fluido (Lukasievicz 2012): tabelas de Talbot contra os limites conhecidos
mf = MODELS["tmfluido"]
tf = np.linspace(0, 0.3, 120)
pf = mf.defaults()
pf.update(pf, thf=0.0)
ham = MODELS["tmhamtl"]
ph = ham.defaults()
for k in ("m", "w0", "v"):
    ph[k] = pf[k]
ph["d"], ph["th"] = pf["d"], pf["th"]
pf0 = dict(pf, cf=1e-9, thf=0.0)   # κ→0: fluido sem condutividade → só a superfície livre
e_ham = np.abs(mf.intensity(tf, pf0) - ham.intensity(tf, ph)).max() / np.abs(ham.intensity(tf, ph) - 1).max()
print(f"TM fluido κ→0 contra TM HAM: erro relativo ao sinal {e_ham:.1e}")
ok &= e_ham < 1e-4

txt = (ROOT / "tlfit" / "examples" / "tmfluido_aco_agua_sintetico.dat").read_text(encoding="utf-8")
dd, _ = parse_table(txt)
pp = mf.defaults()
pp.update(mf.presets["Água"] if "Água" in mf.presets else {})
pp.update(m=32.8, w0=7.4e-5, v=6.9)
true_f = dict(d=6.94e-6, th=-0.062, thf=-26.9)
prob = FitProblem(mf, dd[:, 0], dd[:, 1], pp, [make_free("d", 1e-7, 1e-4, True), make_free("th", -5, 5), make_free("thf", -1e3, 1e3)])
r = run_fit(prob, method="auto", log=lambda s: None)
zs = {k: (r["params"][k] - true_f[k]) / r["sigma"][k] for k in true_f}
print("ajuste   tmfluido  " + "  ".join(f"{k} = {r['params'][k]:.5g} ({zs[k]:+.1f}σ)" for k in true_f)
      + f"  ({r['seconds']:.1f} s)")
ok &= all(abs(z) < 4 for z in zs.values())

# 14. TL + lente de população (Silva 2013): forma fechada em z contra quadratura numérica
from scipy.integrate import quad  # noqa: E402
from tlfit.models import _pl_unit_phase  # noqa: E402

gg = np.array([0.0, 0.05, 0.3, 1.0, 3.0])
tt = np.array([[1e-4], [1e-3], [5e-3], [3e-2]])
worst = 0.0
for ae, l in [(240.0, 3.79e-3), (2000.0, 3.79e-3), (5.0, 3.79e-3)]:   # 5 → ramo de Gauss–Legendre
    h, tau = 40.0, 1.7e-3
    f = _pl_unit_phase(tt, np.exp(-2 * 28.7 * gg), ae, l, h, tau)
    for i, t1 in enumerate(tt[:, 0]):
        for j, u in enumerate(np.exp(-2 * 28.7 * gg)):
            ref = quad(lambda z: u * np.exp(-ae * z) / (1 / tau + h * u * np.exp(-ae * z))
                       * (1 - np.exp(-t1 * (1 / tau + h * u * np.exp(-ae * z)))), 0, l,
                       epsabs=1e-14, epsrel=1e-12)[0]
            worst = max(worst, abs(f[i, j] - ref) / np.abs(f).max())
print(f"TL população: Φ_PL fechada contra quadratura em z: erro relativo máx {worst:.1e}")
ok &= worst < 1e-10

mp = MODELS["tlpop"]
pq = mp.defaults()
tq = np.linspace(0, 0.05, 100)
sh = dict(d=pq["d"], th=0.01, m=pq["m"], w0=pq["w0"], v=pq["v"])
e_sh = np.abs(mp.intensity(tq, dict(pq, th=0.01, tpl=0.0)) - MODELS["tlshen"].intensity(tq, sh)).max()
print(f"θ_PL = 0 contra TL de Shen: diferença máx {e_sh:.1e}")
ok &= e_sh < 1e-6

txt = (ROOT / "tlfit" / "examples" / "tlpop_zblan_sintetico.dat").read_text(encoding="utf-8")
dd, _ = parse_table(txt)
pp = mp.defaults()
true_p = dict(d=2.6e-7, th=0.15, tpl=-1.2e4, tau=1.7e-3)
q = dict(pp, d=1e-7, th=0.05, tpl=-3000.0, tau=5e-4)
prob = FitProblem(mp, dd[:, 0], dd[:, 1], q, [make_free("d", 1e-8, 1e-5, True), make_free("th", -1, 1),
                                               make_free("tpl", -1e5, 1e5), make_free("tau", 1e-4, 1e-1, True)])
r = run_fit(prob, method="auto", log=lambda s: None)
zs = {k: (r["params"][k] - true_p[k]) / r["sigma"][k] for k in true_p}
print("ajuste   tlpop     " + "  ".join(f"{k} = {r['params'][k]:.5g} ({zs[k]:+.1f}σ)" for k in true_p)
      + f"  ({r['seconds']:.1f} s)")
ok &= all(abs(z) < 4 for z in zs.values())

print("OK" if ok else "FALHOU")
sys.exit(0 if ok else 1)
