"""Validação do motor novo.

1. Compara com curvas geradas pelos programas originais (tllam.c e
   tmhamtl-onoff.c, compilados com GSL) — referências em tests/ref/.
   O original imprime 6 casas decimais e integra com tolerância relativa
   1e-4/1e-5, então diferenças de até ~1e-5 são esperadas.
2. Autoconvergência: resultado com quadraturas 2× mais finas.
3. Se c/libtlcore.so (ou tlcore.dll) existir, compara C e Python.

Uso: python tests/validate.py
"""
import ctypes
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
ORDER = {"tllam": "d ths thf m w0 v l df cs cf dens denf ns nf E coefexp dndts nu qpen qpar".split(),
         "tmhamtl": "d th m w0 v beta nu dndt df of".split()}


def load_c():
    for name in ("libtlcore.so", "libtlcore.dylib", "tlcore.dll"):
        f = ROOT / "c" / name
        if f.exists():
            lib = ctypes.CDLL(str(f))
            D = ctypes.POINTER(ctypes.c_double)
            lib.tl_intensity.argtypes = [ctypes.c_int, D, D, ctypes.c_int, D, ctypes.c_double]
            return lib, D
    return None, None


lib, D = load_c()
ok = True
for key, ref, upd in CASES:
    m = MODELS[key]
    p = m.defaults()
    p.update(upd)
    r = np.loadtxt(REF / ref)
    t = np.ascontiguousarray(r[:, 0])
    I = m.intensity(t, p)
    I2 = m.intensity(t, p, Accuracy().scaled(2))
    d_orig, d_self = np.abs(I - r[:, 1]).max(), np.abs(I - I2).max()
    line = f"{key:8s} {ref:34s} |novo−original| {d_orig:.1e}   autoconvergência {d_self:.1e}"
    if lib:
        pa = np.array([p[k] for k in ORDER[key]], float)
        out = np.empty_like(t)
        lib.tl_intensity(0 if key == "tllam" else 1, pa.ctypes.data_as(D), t.ctypes.data_as(D),
                         len(t), out.ctypes.data_as(D), 1.0)
        line += f"   |C−Python| {np.abs(out - I).max():.1e}"
    print(line)
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

dat = np.loadtxt(ROOT / "tlfit" / "examples" / "tllam_amostra_agua.dat")
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

print("OK" if ok else "FALHOU")
sys.exit(0 if ok else 1)
