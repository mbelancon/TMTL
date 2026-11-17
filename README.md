# TMTL FIT - Fast fitting

*English · [Português](README.pt-BR.md)*

Fitting of time-resolved thermal lens and thermal mirror transients, with a graphical
interface and fitting of transient series (several powers). The interface is in English and
Portuguese and follows the system language (there is a PT/EN switch at the top).

> **Portable version for Windows, no Python needed:** download the zip
> `TMTL2-…-portable-windows.zip` from the
> [Releases](https://github.com/mbelancon/TMTL/releases) page, unzip it and double-click
> `Open TMTL2.bat`. Details in [Usage](#usage).

## Models

Names: **TM** is the thermal mirror and **TL** the thermal lens. When the model adds the TL or
the TM of a solid to the thermal lens of the fluid in contact with it (air, water), the name gets
"+TL-Fluid": **TL+TL-Fluid** and **TM+TL-Fluid**. The fluid always contributes a thermal lens;
there is no thermal mirror in the fluid.

The acronyms LAM, BAM and HAM (Sato et al., 2008) stand for *Low*, *Beer* and *High absorption
model*: **LAM** = Low Absorption Model (A_e·thickness ≪ 1), **BAM** = Beer Absorption Model
(absorption by Beer's law, valid for any A_e) and **HAM** = High Absorption Model (heat deposited
at the surface).

| Model | Use | Reference |
|---|---|---|
| **TL+TL-Fluid (TLLAM)** | thermal lens of the sample plus that of the fluid, with sample–fluid heat coupling and any thickness (unified model) | Malacarne et al., *Appl. Spectrosc.* **65**, 99 (2011) |
| **Shen TL (classic)** | classic closed-form expression of the thermal lens in the mode-mismatched configuration | Baesso, Shen & Snook, *J. Appl. Phys.* **75**, 3732 (1994), eq. 10 |
| **TL with photobleaching** | thermal lens of a liquid whose absorption decays by a photochemical reaction (rate K_T) down to ε̄·β; fits θ, K_T and ε̄ | Herculano et al., *J. Phys. Chem. B* **117**, 1932 (2013), eqs. 3–6 |
| **TL + population lens** | thermal lens plus excited-state lens (saturable population, lifetime τ); fits θ_th, θ_PL, τ and D | Silva et al., *Appl. Phys. Lett.* **102**, 141910 (2013), eqs. 1–3 |
| **TM LAM** | thermal mirror, low absorption (source uniform in z), semi-infinite sample | Sato et al. (2008), eq. 32 |
| **TM LAM, finite sample** | low-absorption thermal mirror in a sample of thickness *l* (free faces) | Astrath et al., *J. Opt. Soc. Am. B* **28**, 1735 (2011), eq. 10 |
| **TM BAM** (Beer's law) | thermal mirror for any A_e; tends to LAM and to HAM | Sato et al. (2008), eqs. 33–34 |
| **TM+TL-Fluid, HAM on/off (approximation: no fluid feedback)** | thermal mirror in an opaque sample, with switch-off, plus the thermal lens of the fluid (air); the sample surface is free, as in `tmhamtl-onoff.c` | Sato et al., *J. Appl. Phys.* **104**, 053520 (2008); Flizikowski et al., *Opt. Express* **28**, 7116 (2020) |
| **TM+TL-Fluid, HAM (full coupling)** | thermal mirror with the thermal lens of the fluid in contact (air, water); sample–fluid heat coupling by Laplace transform (Talbot tables); fits D, θ and θ_f | Lukasievicz et al., *Appl. Spectrosc.* **66**, 1461 (2012); Flizikowski et al., *Opt. Express* **28**, 7116 (2020), eqs. 3–10 |

**Which TM+TL-Fluid to use.** Both include the heat that flows from the sample into the fluid
and the thermal lens it forms. The difference is the feedback: in **HAM on/off (approximation)**,
equivalent to `tmhamtl-onoff.c`, the sample surface is free, and the fluid receives heat but does
not cool the sample back; in **HAM (full coupling)** (eqs. 3–10 of Flizikowski et al., 2020) the
fluid also cools the sample. Without that feedback the two agree to 10⁻¹¹. For the three metals of
the paper (copper, bronze and stainless steel, in air) the approximation differs from the full
model by 10⁻⁴ to 10⁻³ of the signal with the laser on and by 10⁻³ to 10⁻² right after switch-off
(1.3 % for steel). To fit data from the paper, or when precision right after switch-off matters,
use the full coupling; to compare with the original program, the approximation.

**Each model has its own θ.** The interface shows the definition and unit of the chosen model's θ
right below its description, and the report and the series repeat that definition:

| Model | Fitted θ | Unit |
|---|---|---|
| TL+TL-Fluid (TLLAM) | sample θ = ϑ·10⁻⁶, ϑ = P_e A_e φ/(k_s λ_p) (the phase is ϑ·l·χ(α,l)) | 10⁶ K/m |
| TL+TL-Fluid (TLLAM) | fluid θ = 2 P_e A_e φ (∂n_f/∂T)/λ_p | W/(K·m²) |
| TM+TL-Fluid, HAM on/off (approximation) | θ_HAM = θ_TM/A_e = −P_e α_T (1+ν) φ/(λ_p k) | dimensionless |
| TM LAM, TM LAM finite, TM BAM | θ_TM = −P_e A_e α_T (1+ν) φ/(λ_p k) | 1/m |
| Shen TL | θ = −P_e A_e L φ (ds/dT)/(k λ_p) | rad |
| TL with photobleaching | θ = −P_e β L φ (dn/dT)/(k λ_p), β = initial absorption (the final plateau is ε̄·θ) | rad |
| TM+TL-Fluid, HAM (full coupling) | θ = θ_HAM of the sample (dimensionless); θ_f = 2 P_e φ (dn/dT)_f/λ_p, in W/(K·m) | — |
| TL + population lens | θ_th = −P_abs (ds/dT) φ/(k λ_p), P_abs = P A_e L; θ_PL = I0 σ_abs N_T Δα_P (2π f_L)²/(n λ_p E) | rad; 1/(m·s) |

The thermal mirror models accept switching off the excitation at `t_off` (0 = no switch-off). The
sign convention of the phase is the same as in `tmhamtl-onoff.c` (phase added to V·g in the
diffraction integral); relative to eq. 31 of Sato et al. this amounts to flipping the sign of
θ_TM, with no other effect on the fit.

**Relation to θ_TM of Flizikowski et al. (2020).** In the paper, θ_TM = −P_e (1−R) α_T (1+ν) φ/(λ_p k),
with P_e the incident power and R the reflectivity. It is the θ_HAM of the table above with the
absorbed power P_e(1−R) in place of P_e; that is, with incident P_e, θ_TM(paper) = (1−R)·θ_HAM.
Table 1 of the paper lists the magnitude of θ_TM/P_e (by the formula it is negative for a sample
that expands, as here).

In the photobleaching model, fix D with a measurement of the pure solvent (K_T = 0, ε̄ = 1) and fit
θ, K_T and ε̄ on the transients of the solution. In a series (several P_e), K_T = K_m + K_pb·P_e:
the line of K_T against power gives K_m (mass diffusion) and K_pb (photobleaching).

The Shen model is first order in θ: for |θ| ≳ 0.1 rad the difference from the full diffraction
integral exceeds ~1 % of the signal. For large signals, use the integral models.

The examples (`tlfit/examples/`) carry a line `# tlfit: model=… m=… w0=…` that, when the example
is opened in the interface, selects the model and fills in the geometry (the older `modelo=` is
also accepted).

## Usage

### Portable version for Windows: no Python needed

On the repository's **Releases** page, download `TMTL2-<version>-portable-windows.zip`, unzip it
wherever you like and double-click **`Open TMTL2.bat`**. The folder already contains the
official Python (the "embeddable" version from python.org), numpy and scipy: **nothing needs to be
installed**, no administrator rights are needed, and nothing on the computer is changed.

The interface opens in the browser, in the system language. To quit, close the tab (the program
exits by itself about 30 s later) or the black window. `Check installation.bat` checks the
calculations.

### From the source code: if you already have Python

```bash
pip install numpy scipy          # Anaconda/Miniconda already has them
python app.py                    # or double-click rodar.bat (finds Anaconda by itself)
```

Or install it as a package, straight from GitHub, and open it with `tlfit`:

```bash
pip install git+https://github.com/mbelancon/TMTL
tlfit
```

The interface language comes from the browser, and the choice made with the PT/EN switch is
remembered.

## Transient series (several powers)

In the **Series** selector, at the top of the sidebar:

1. **Open folder** (or drag the folder, or **Choose files**). The files `.dat`, `.txt`, `.csv`
   and `.tsv` are read, using the first two columns (time and signal).
2. The power of each file is read from its name (`glass_150mW.dat`, `sample 0,2W.txt`,
   `P150.dat`) or from the header (`# Power: 150 mW`, `Potência = 0.15 W`) and converted to mW.
   The table shows the value and lets you correct it or exclude files. The name of the variable
   can be changed (concentration, temperature, oxidation time, etc.).
3. If the signals are not normalized (in volts, for example), normalization by I(t=0) is turned
   on automatically; it is worth leaving **I0 (scale)** free in that case.
4. **Fit series** fits all transients with the same model, the same fixed parameters and the
   same method (in automatic mode, ~2 s per transient with TLLAM).

Results:

* transients and fits overlaid, coloured by power (click a table row to highlight a transient);
* any fitted parameter × power, with error bars (1σ):
  * for the θ (proportional to power): a line through the origin, with slope θ/P_e, and a line
    with intercept, which warns when the intercept differs from zero; the "θ / power" option shows
    θ/P of each transient and the mean;
  * for D, I0 and the others: mean, standard deviation, standard error and the slope against
    power, with a warning when it is significant (heating, non-linear effects);
* **Export results (.csv)**: one row per transient with values, uncertainties, θ/P, χ² and R²;
  **Series report (.txt)**: configuration, results and statistics.

The folders `tlfit/examples/serie_vidro_ar` (normalized) and `serie_vidro_volts` (in volts)
contain synthetic series for testing.

## Building the portable version

* **Automatic, on GitHub:** the workflow `.github/workflows/build.yml` builds and tests the
  portable version on a GitHub Windows runner at every push (the zip is under *Actions → run →
  Artifacts*). When a tag is created, for example `git tag v2.2.1 && git push --tags`, the zip is
  published under **Releases**.
* **On your computer (Windows, with internet):** double-click `montar_portatil.bat`. The script
  `ferramentas/montar_portatil.ps1` downloads the embeddable Python from python.org, installs
  numpy and scipy inside it (using the pip of a Python you already have, only to download the
  packages), copies the program, runs the numerical validation with the embedded Python and
  builds the zip in `%USERPROFILE%\TMTL2-portable` (outside synchronized folders: it is ~150 MB).
  The log goes to `montar_portatil.log`.

The versions of Python, numpy and scipy are pinned in the script, so the portable version is
reproducible.

## Shared server (optional)

The interface can also be hosted for several people: `python -m tlfit --host 0.0.0.0
--no-browser` on a machine in the network, or the `Dockerfile` (port 7860) on a service such as
Hugging Face Spaces. `TLFIT_MAX_FITS` (default 2) limits simultaneous fits; the others wait in the
queue. There is no login. Each user sees the interface and the messages in their own language.

## Fitting method

**Automatic (default): two-level fit.** The θ and the scale I0 enter the phase linearly,
Φ = θ_sample·Φ_s + θ_fluid·Φ_f (in the TM, Φ = θ·Φ₁). Therefore:

1. the other free parameters (typically only D; at most two) are swept on a grid: 40 points on a
   log scale for one parameter, 12 × 12 for two;
2. at each point, Φ_s and Φ_f are computed only once, and the θ (and I0) are fitted inside,
   starting from several points, at a cost 10–15× lower than a full model evaluation, with
   analytic derivatives;
3. the best local minima of the grid are refined by Levenberg–Marquardt with all parameters
   free, at normal accuracy, which also gives the uncertainties.

Because it sweeps the whole range of the non-linear parameters, the result does not depend on
luck, and the log reports how many local minima the sweep found. The "Search" panel shows the
points of the sweep (for example, optimal θ × D, with the colour giving χ²). With more than two
free parameters besides the θ and I0, the automatic mode uses the swarm:

* **PSO** (particle swarm, inertia 0.9 → 0.4, global topology) for the global search, with the
  user's initial guess as one of the particles and the others in a Latin hypercube. Positive
  parameters whose range spans more than two decades (diffusivities) are explored on a log scale.
  The PSO + LM pair is repeated in **runs** (default: 3) with different seeds, and the best is
  kept; the log warns if the runs disagree.
* **Levenberg–Marquardt** (trust region with bounds, `scipy.optimize.least_squares`) from the best
  point. It converges in a few evaluations and provides the Jacobian, from which come the
  uncertainties (1σ, scaled by the reduced χ²) and the correlation matrix. The log warns when a
  parameter ends at a bound or when two parameters have a correlation above 0.95 (for example,
  sample θ × fluid θ in the air example).

The search (sweep or swarm) uses a lighter numerical mesh: g up to 18 instead of 30 and panels
0.7×, with an error in I(t) below 3e-7, against a typical noise of ~1e-4. The final refinement
always uses normal accuracy.

Times of a full fit (2 cores):

| Example | PSO + LM, 3 runs | Automatic |
|---|---|---|
| TLLAM, sample in air | 4.7 s | 2.0 s |
| TLLAM, sample in water | 4.9 s | 1.8 s |
| TM on/off | 22 s | 5.9 s |
| TM on/off + lens in air, I0 free | ~50 s (estimated) | 16 s |

Both modes reach the same parameters in these examples.

What the original notebook did was a random search around the best point (perturbation of
±percent %), not a PSO: there was no velocity or particle memory, and the reach of the search was
limited by the initial guess.

## Mathematical methods

The document [`docs/MATH_METHODS.md`](docs/MATH_METHODS.md) (in Portuguese:
[`docs/METODOS_MATEMATICOS.md`](docs/METODOS_MATEMATICOS.md)) describes the general problem
(heating, phase of the probe beam and diffraction up to the detector), shows that all models share
the same diffraction integral and differ only in the phase kernel, gives the equations of each
model and explains the techniques used to speed up the integrations.

## What changed in the computation

The model is

```math
I(t) = (1+V^2)\left|\int_0^{30} e^{-g}\, e^{\,i\,[V g + \Phi(g,t)]}\, dg\right|^2,\qquad
\Phi(g,t) = \int_0^\infty K(\alpha,t)\left[J_0\!\left(\alpha\, w_0\sqrt{m g}\right) - 1\right] d\alpha
```

and in the original code each Φ was an adaptive integral (QAG) called **inside** another adaptive
integral, twice (once for cos, once for sin), and with the lens in air there was a third level
(an integral in τ) inside the integral in α — also redone for each g, although it does not depend
on g.

The changes:

1. **Separation of variables.** The Bessel factor does not depend on t or on the thermal
   parameters. With fixed quadratures in α and g, Φ becomes a matrix product
   `Φ = K(t,α) · M(α,g)`, and `M` is computed once per geometry (m, w₀).
2. **A single complex sum** gives the real and imaginary parts together.
3. **Integral in τ outside the loop in g**, with the substitution τ = t sin²φ (removes the
   τ^{-1/2} singularity) and tanh-sinh quadrature.
4. **Panelled Gauss–Legendre quadratures**, with geometric panels near α = 0 and automatic
   refinement in g where the total phase turns fast (large θ). Self-convergence (2× resolution)
   is ~1e-9 or better.
5. Stable forms: `(cosh x − 1)/(x + sinh x)` without overflow, `expm1` for `1 − e^{−x}`.

## Validation

`python tests/validate.py` compares with curves generated by the original programs (in
`tests/ref/`), checks self-convergence, compares the fluid term with the numerical inversion of
the Laplace solution in both regimes of the formula (erf and Dawson), fits the synthetic examples
and reproduces the paper by Flizikowski et al. (2020):

    tllam    tllam_d5.317e-7_ths4.33_thf0.dat   |new−original| 4.9e-07   self-convergence 2.3e-09
    tllam    tllam_d5e-7_ths4_thf30.dat         |new−original| 5.0e-07   self-convergence 2.8e-09
    tmhamtl  tmhamtl_default.dat                |new−original| 1.5e-04   self-convergence 4.5e-14
    fluid    Laplace, erf branch: max relative error 2.3e-08
    fluid    Laplace, Dawson branch: max relative error 3.6e-11
    fit      water, automatic mode: χ² = 9.05310e-06, D = 4.97742e-07, θs = 4.6501, θf = -1197.1 (2.2 s)
    series   θs/P = 0.02005 ± 0.00021 rad/mW (true 0.020); D = 5.3402e-07 ± 1.2e-08 (true 5.4e-7)
    TM-BAM   Ĝ(α,t) against quadrature in τ: max relative error 3.4e-08
    TM-LAM   phase against eq. 32 of Sato et al. (2008): max relative error 2.4e-15
    limits   BAM→HAM 1.7e-09   BAM→LAM 2.1e-08   finite→LAM 1.3e-10 (relative to the signal)
    TL Shen  against the diffraction integral (θ = 0.01): difference 1.1e-03 of the signal (order θ)
    fit      tmlam     D = 6.0163e-07 (+1.1σ)  θ = -120.02 (-0.3σ)  (4.2 s)
    fit      tmbam     D = 6.0266e-07 (+2.0σ)  θ = -149.9 (+2.1σ)  (7.3 s)
    fit      tmfinito  D = 5.0145e-07 (+1.4σ)  θ = -149.93 (+1.0σ)  (4.1 s)
    fit      tlshen    D = 4.9652e-07 (-1.3σ)  θ = -0.10419 (-1.3σ)  (0.7 s)
    TL bleach phase against eq. 4 of Herculano et al. (2013): max relative error 8.1e-12
    limits   K_T→0 8.8e-04   ε̄=1 8.8e-04 (relative to the signal; Shen is of order θ)
    fit      tlbranq   th = 0.048085 (+1.8σ)  kt = 2.0064 (+1.1σ)  eb = 0.27033 (+0.8σ)  (12.3 s)
    TM fluid κ→0 against TM HAM: error relative to the signal 7.6e-12
    fit      tmfluido  d = 6.9307e-06 (-0.5σ)  th = -0.062098 (-0.7σ)  thf = -26.858 (+0.5σ)  (19.2 s)
    TL population: closed Φ_PL against quadrature in z: max relative error 4.2e-15
    θ_PL = 0 against Shen TL: max difference 7.8e-07
    fit      tlpop     d = 2.5992e-07 (-0.1σ)  th = 0.14999 (-0.1σ)  tpl = -12071 (-0.5σ)  tau = 0.0016903 (-0.4σ)  (76.4 s)
    paper Flizikowski 2020, copper : TM+TL-Fluid on 2.4e-05, off 2.4e-05   (TM+TL-Fluid HAM on/off, approximation without fluid feedback: 9.2e-04; relative to the signal)
    paper Flizikowski 2020, bronze : TM+TL-Fluid on 3.5e-06, off 5.2e-06   (TM+TL-Fluid HAM on/off, approximation without fluid feedback: 4.3e-03; relative to the signal)
    paper Flizikowski 2020, steel  : TM+TL-Fluid on 1.5e-06, off 1.3e-06   (TM+TL-Fluid HAM on/off, approximation without fluid feedback: 1.3e-02; relative to the signal)
    i18n     model texts without a Portuguese translation: 0
    OK

The reference cases from the original programs are in air, where the original code is correct.
The times in parentheses depend on the computer. The last check makes sure every model text shown
in the interface has a Portuguese translation in `tlfit/i18n.py`.

Eqs. 3–10 of Flizikowski et al. (2020) (G of eq. 4, u_z of eq. 6, lens in air of eq. 8, signal of
eq. 10) were implemented separately, with the parameters of Table 1 and fig. 1 of the paper, and
TM+TL-Fluid (full coupling) reproduces them. Switch-off is the superposition T_on(t) − T_on(t − ξ).
`tmhamtl` (free surface, as in `tmhamtl-onoff.c`) differs from the paper by 9e-4 (copper), 4e-3
(bronze) and 1.3e-2 (steel) of the signal in the off regime, because the paper includes the
sample–air heat exchange; without that exchange (κ → 0) the two models agree to 1e-11. The
difference grows with k_air/k_sample.

* The BAM kernel, Ĝ(α,t) = ∫T̂(α,z,t)e^{−αz}dz, comes from
  dĜ/dt = [A_e erfc(αs) − α e^{−α²s²} erfcx(A_e s)]/(A_e² − α²), s = √(Dt), integrated in closed
  form. Where the closed form loses precision (small s·max(α, A_e), or α ≈ A_e) a series in s is
  used in which the factor (A_e² − α²) cancels term by term. The result is compared with direct
  numerical quadrature of dĜ/dt at 400 random points (0.1 < A_e < 10⁸ m⁻¹, 10⁻⁶ < t < 10 s,
  including α ≈ A_e).
* LAM agrees with the closed form in g of eq. 32 of Sato et al.; BAM tends to the HAM of
  `tmhamtl-onoff.c` (including switch-off) when A_e grows, and to LAM when A_e → 0.
* Fitting curves generated by BAM with LAM (geometry of fig. 11 of Sato et al.), the errors in
  θ_TM and t_c are ~1.5 % for A_e = 77 m⁻¹ and 16–18 % for 1061 m⁻¹, as in the paper: above
  ~200 m⁻¹ use BAM.
* The fits of the synthetic examples recover D and θ within ±2.1σ.

For TLLAM the difference is the rounding of the original's output (6 decimal places). For the TM
the difference comes from the original's tolerance (QAG with relative error 1e-4): recompiling the
original program with tolerance 1e-8, the points of largest discrepancy (t = 0.37 s and 1.85 s)
agree with the new code in all 6 printed decimal places (0.790799 and 0.998527).

## Use as a library

```python
import numpy as np
from tlfit.models import MODELS
from tlfit.optim import FitProblem, make_free, run_fit

t, y = np.loadtxt("datatofit.dat", unpack=True)
m = MODELS["tllam"]
p = m.defaults(); p.update(d=5e-7, ths=4, thf=0)
prob = FitProblem(m, t, y, p, [make_free("d", 1e-8, 1e-5), make_free("ths", 0, 20)])
res = run_fit(prob)            # automatic mode (two levels)
print(res["params"]["d"], res["sigma"]["d"])
```

Log messages are `tlfit.i18n.Msg` objects (they behave as English strings);
`tlfit.i18n.tr(msg, "pt")` gives the Portuguese text.

## Structure

    app.py              entry point (opens the interface)
    rodar.bat           runs from the source code (finds Anaconda/Miniconda)
    montar_portatil.bat builds the portable version (ferramentas/montar_portatil.ps1)
    .github/workflows/  automatic build of the portable version and Releases
    legado/             original code (C with GSL), notebook and README of version 1
    tlfit/models.py     models (TL+TL-Fluid, Shen TL, TL with photobleaching, TL + population, TM LAM, TM LAM finite, TM BAM, TM+TL-Fluid HAM on/off and full coupling)
    tlfit/optim.py      two-level fit, PSO, Levenberg–Marquardt, uncertainties
    tlfit/batch.py      series fitting and statistics (θ/P, means, trends)
    tlfit/server.py     local server (standard library) for the interface
    tlfit/i18n.py       Portuguese translations of the messages generated in Python
    tlfit/static/       interface (HTML/JS, no external dependencies, in English and Portuguese)
    tlfit/examples/     example transients and series
    tests/              validation against the original code
    docs/               mathematical methods (model equations and numerical techniques), in English and Portuguese
    CLAUDE.md           rules for working on the code with Claude Code (imports CONTEXTO.md; in Portuguese)
    CONTEXTO.md         short context for whoever works on the code (decisions, conventions; in Portuguese)
    CHANGELOG.md        changes per version

## References

Papers used or cited in the project, in chronological order. The Models table says which equation
of each one was implemented.

1. Shen, Lowe & Snook, *Chem. Phys.* **165**, 385 (1992). Mode-mismatched thermal lens model, the
   basis of the Shen TL.
2. Baesso, Shen & Snook, *J. Appl. Phys.* **75**, 3732 (1994). Shen TL (eq. 10).
3. Sato et al., *J. Appl. Phys.* **104**, 053520 (2008). TM LAM (eq. 32), TM BAM (eqs. 33–34) and
   the HAM limit of TM+TL-Fluid, HAM on/off.
4. L. C. Malacarne, N. G. C. Astrath, G. V. B. Lukasievicz, E. K. Lenzi, M. L. Baesso and
   S. E. Bialkowski, "Time-Resolved Thermal Lens and Thermal Mirror Spectroscopies with
   Sample-Fluid Heat Coupling: A Complete Model for Material Characterization",
   *Appl. Spectrosc.* **65**, 99 (2011), doi:10.1366/10-06096. TL+TL-Fluid (TLLAM).
5. N. G. C. Astrath, L. C. Malacarne, V. S. Zanuto, M. P. Belançon, R. S. Mendes, M. L. Baesso and
   C. Jacinto, "Finite-size effect on the surface deformation thermal mirror method",
   *J. Opt. Soc. Am. B* **28**, 1735 (2011). TM LAM, finite sample (eq. 10).
6. G. V. B. Lukasievicz, L. C. Malacarne, N. G. C. Astrath, V. S. Zanuto, L. S. Herculano and
   S. E. Bialkowski, "A Theoretical and Experimental Study of Time-Resolved Thermal Mirror with
   Non-Absorbing Heat-Coupling Fluids", *Appl. Spectrosc.* **66**, 1461–1467 (2012).
   TM+TL-Fluid, HAM, full coupling (eqs. 4–9 and 18–22).
7. Herculano et al., *J. Phys. Chem. B* **117**, 1932 (2013). TL with photobleaching (eqs. 3–6).
8. Silva et al., *Appl. Phys. Lett.* **102**, 141910 (2013). TL + population lens (eqs. 1–3).
9. G. A. S. Flizikowski, O. A. Capeloto, V. G. Camargo, B. Anghinoni, M. L. Baesso,
   L. C. Malacarne, M. P. Belançon, T. Požar and N. G. C. Astrath, "Laser induced thermoelastic
   surface displacement in solids detected simultaneously by photothermal mirror and
   interferometry", *Opt. Express* **28**, 7116 (2020), doi:10.1364/OE.386344. Thermal mirror with
   on/off excitation and thermal lens in air, with sample–fluid heat coupling (eqs. 3–10), in
   copper, bronze and stainless steel. Reproduced in `tests/validate.py` (see Validation). The
   default values of *m* and *w₀* of TM+TL-Fluid, HAM on/off match the geometry of the paper
   (w_e = 646 µm, w_p = 1584 µm).

## License

Copyright © Marcos P. Belançon and contributors.

This program is free software, distributed under the **GNU General Public License, version 3 or
later** (GPL-3.0-or-later): it can be used, studied, modified and redistributed, provided that
modified versions are also distributed under the GPL. The full text is in [`LICENSE`](LICENSE).

The portable version includes Python (PSF license), numpy and scipy (BSD licenses); the texts of
those licenses come with each component inside the `python` folder of the portable version.

## Acknowledgements

I thank the State University of Maringá (Universidade Estadual de Maringá), especially the
Graduate Program in Physics. Special thanks to the professors who led the projects from which the
models implemented here resulted: Luis C. Malacarne
([Google Scholar](https://scholar.google.com.br/citations?user=AeaidPoAAAAJ&hl=pt-BR)),
Nelson G. C. Astrath
([Google Scholar](https://scholar.google.com.br/citations?user=9ryRh9oAAAAJ&hl=pt-BR)) and
Mauro L. Baesso
([Google Scholar](https://scholar.google.com.br/citations?user=Kea32YMAAAAJ&hl=pt-BR)).
