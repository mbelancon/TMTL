# Mathematical methods

*English · [Português](METODOS_MATEMATICOS.md)*

This document describes what the `tlfit` models compute, why the **phase** is the only thing that
changes from one model to another, and which numerical techniques make the computation fast. The
equations follow the code in [`tlfit/models.py`](../tlfit/models.py) (the notation in backticks is
the code's); the references of each model are in the Models table of the
[README](../README.md#models).

Naming convention: **TL** = thermal lens, **TM** = thermal mirror, **+TL-Fluid** = thermal lens of
the fluid in contact with the sample (the fluid only contributes a TL; there is no TM in the
fluid). LAM/BAM/HAM = *Low/Beer/High Absorption Model*.

Contents

1. [The general problem](#1-the-general-problem)
2. [The phase is what differentiates the models](#2-the-phase-is-what-differentiates-the-models)
3. [The particularity of each model](#3-the-particularity-of-each-model)
4. [Techniques that speed up the integrations](#4-techniques-that-speed-up-the-integrations)
5. [Techniques in the fit](#5-techniques-in-the-fit)

---

## 1. The general problem

All models describe the same experiment, in three stages:

1. **Heating.** The excitation beam, Gaussian with radius $`w_{0e}`$ at the sample (`w0`), deposits
   heat with radial profile $`e^{-2r^2/w_{0e}^2}`$ and a depth profile $`Q(z)`$ that depends on the
   absorption. The temperature $`T(r,z,t)`$ obeys the diffusion equation

   ```math
   \frac{\partial T}{\partial t} - D\,\nabla^2 T = Q_0\, e^{-2r^2/w_{0e}^2}\, Q(z),
   ```

   with $`D`$ the thermal diffusivity (`d`) and the boundary conditions of each model (free
   surface, finite sample, heat exchange with the fluid).

2. **Phase change of the probe beam.** Heating changes the optical path: the refractive index
   changes with temperature (thermal lens, in the sample or in the fluid) and the surface deforms
   by expansion (thermal mirror). The probe beam, of radius $`w_{1p}`$ at the sample, leaves with a
   phase $`\Phi(r,t)`$ measured relative to the centre, $`\Phi(r,t) = \phi(r,t) - \phi(0,t)`$.

3. **Propagation to the detector.** The probe beam propagates (Fresnel diffraction) to the
   detector, in the far field, which measures only the centre of the beam. With the radial
   variable $`g = (r/w_{1p})^2`$, the normalized intensity is

   ```math
   I(t) = (1+V^2)\,\left|\int_0^{\infty} e^{-g}\, e^{\,i\,[V g + \Phi(g,t)]}\, dg\right|^2 ,
   ```

   where $`m = (w_{1p}/w_{0e})^2`$ and $`V`$ (`m`, `v`) are the geometric parameters of the
   mode-mismatched setup ($`V \approx Z_1/Z_c`$, with $`Z_1`$ the distance from the probe waist to the
   sample and $`Z_c`$ the confocal length). Without phase, $`\int_0^\infty e^{-(1-iV)g}dg = 1/(1-iV)`$
   and the factor $`(1+V^2)`$ makes $`I = 1`$: the signal is normalized by its value before heating.
   The code integrates $`g`$ up to $`g_\text{max} = 30`$ (the neglected term is of order $`e^{-30}`$),
   like the original program.

The scale `scale` (I0) multiplies $`I(t)`$ and serves signals that are not exactly normalized.

## 2. The phase is what differentiates the models

Stages 1 and 3 are the same in all models (apart from $`Q(z)`$ and the boundaries). What changes is
**how temperature becomes phase**:

| Effect | Phase (relative to the centre) |
|---|---|
| **TL** (sample of thickness $`L`$) | $`\Phi = \dfrac{2\pi}{\lambda_p}\dfrac{ds}{dT}\displaystyle\int_0^L [T(r,z,t)-T(0,z,t)]\,dz`$ |
| **TM** (reflecting surface) | $`\Phi = \dfrac{2\pi}{\lambda_p}\, 2\,[u_z(r,0,t)-u_z(0,0,t)]`$, with $`u_z \propto -(1+\nu)\alpha_T \displaystyle\int_0^\infty \hat T(\alpha,z,t)\,e^{-\alpha z}dz`$ |
| **TL-Fluid** | $`\Phi_f = \dfrac{2\pi}{\lambda_p}\dfrac{\partial n_f}{\partial T}\displaystyle\int_\text{fluid} [T_f(r,z,t)-T_f(0,z,t)]\,dz`$ (two passes in the TM) |
| Population lens | index change caused by the excited-state population (not thermal) |

In the TL the phase comes from the temperature **integrated over the thickness**; in the TM it
comes from the **surface displacement**, which weights the temperature with $`e^{-\alpha z}`$
(quasi-static thermoelastic solution). That is why the same temperature profile gives different
phases in TL and TM.

### Common form: the kernel $`K(\alpha,t)`$

Writing the temperature by its Hankel transform, $`T(r,z,t) = \int_0^\infty \hat T(\alpha,z,t)\,J_0(\alpha r)\,\alpha\,d\alpha`$,
and using $`r = w_{1p}\sqrt{g} = w_{0e}\sqrt{m g}`$, all thermal phases take the same form:

```math
\Phi(g,t) = \int_0^\infty K(\alpha,t)\,\left[J_0\!\left(\alpha\, w_{0e}\sqrt{m g}\right) - 1\right] d\alpha .
```

The factor $`[J_0 - 1]`$ is the "geometric" part (where the probe sees the sample) and is the same in
all models. **All the physics of each model is in the kernel $`K(\alpha,t)`$**: the source, the
boundaries, the optical effect (TL, TM, TL-Fluid) and the switch-off. The Hankel transform of the
Gaussian source gives the factor $`e^{-\alpha^2 w_{0e}^2/8}`$ that appears in every kernel. In the
code, each model implements only the method `kernel(t, a, p, acc)`; the diffraction integral
(`Model.intensity`) is shared.

Two practical consequences:

- the parameters $`\theta`$ enter $`K`$ **linearly**, and therefore $`\Phi`$:
  $`\Phi = \theta_\text{sample}\Phi_s + \theta_\text{fluid}\Phi_f`$ (in the TM, $`\Phi = \theta\,\Phi_1`$).
  This is used in the two-level fit (section 5);
- switching off the excitation at $`t_\text{off}`$ (`of`) is the superposition
  $`K_\text{off}(t) = K(t) - K(t-t_\text{off})`$, because the heat equation is linear.

Exceptions: the **Shen TL** is a closed expression for $`I(t)`$ (it does not go through the numerical
integral), and the **population lens** has a phase that does not have the form $`K\cdot[J_0-1]`$,
computed directly in $`g`$ (`direct_phase`).

## 3. The particularity of each model

Notation: $`a=\alpha`$, $`w = w_{0e}`$, $`\gamma(\alpha) = e^{-\alpha^2 w^2/8}`$, $`b = D\alpha^2`$ and
$`h(x) = (\cosh x - 1)/(x + \sinh x)`$ (finite thickness factor, `_hh`).

### Summary

| Model | Source $`Q(z)`$ | Boundaries | Phase | Kernel $`K(\alpha,t)`$ |
|---|---|---|---|---|
| TL+TL-Fluid (TLLAM) | uniform (low absorption) | thickness $`l`$, heat exchange with the fluid | sample TL + fluid TL | closed; fluid with erf, erfcx, Dawson |
| Shen TL | uniform | thin sample, no fluid | TL | — (closed form of $`I(t)`$) |
| TL with photobleaching | uniform, decaying in time | thin sample | TL | closed |
| TL + population lens | Beer | thin sample | TL + population | closed TL; population directly in $`g`$ |
| TM LAM | uniform | semi-infinite, free surface | TM | closed |
| TM LAM, finite sample | uniform | thickness $`l`$, free faces | TM | closed, with $`h(l\alpha)`$ |
| TM BAM | Beer, $`A_e e^{-A_e z}`$ | semi-infinite | TM | closed + series near $`\alpha s\to 0`$ |
| TM+TL-Fluid, HAM on/off (approximation) | at the surface | semi-infinite, free surface | TM + fluid TL (air) | closed TM; fluid by quadrature in $`\tau`$ |
| TM+TL-Fluid, HAM (full coupling) | at the surface | heat exchange with the fluid | TM + fluid TL | Laplace–Hankel inverted numerically |

### 3.1 TL+TL-Fluid (TLLAM) — Malacarne et al. (2011)

Thermal lens of the sample (thickness $`l`$, `l`) plus that of the fluid, with heat coupling. The
kernel is the sum of two terms linear in $`\theta_s`$ (`ths`) and $`\theta_f`$ (`thf`):

```math
K = \theta_s\, l\,\gamma(\alpha)\,B(\alpha)\,\frac{1-e^{-D\alpha^2 t}}{\alpha} \;+\; \theta_f\,K_f(\alpha,t).
```

The factor $`B(\alpha)`$ (the $`\chi`$ of the README) gathers the three optical contributions of the
sample: $`dn/dT`$, the deformation of both faces and the photoelastic effect,

```math
B = \frac{dn}{dT} + 4(n_s-n_f)(1+\nu)\alpha_T\frac{h(l\alpha)}{l\alpha}
  + \frac{n_s^3 E\alpha_T}{4(1-\nu)}\left[q_\parallel + 3q_\perp - 4\big(q_\parallel\nu + q_\perp(2+\nu)\big)\frac{h(l\alpha)}{l\alpha}\right].
```

The fluid term, with $`k = D c_s\rho_s`$, $`k_f = D_f c_f \rho_f`$ and
$`t_0 = (k^2-k_f^2) D D_f/(k^2 D_f - k_f^2 D)`$, is

```math
K_f = \frac{\gamma(\alpha)}{(k^2-k_f^2)\alpha^2}\Big[k\,\mathrm{erf}(\alpha\sqrt{D_f t}) - k_f\,\mathrm{erf}(\alpha\sqrt{D t})
      + \mathcal{T}(D_f, k) - \mathcal{T}(D, k_f)\Big],
```

with, for each pair $`(\delta, c)`$,

```math
\mathcal{T}(\delta,c) =
\begin{cases}
c\sqrt{\dfrac{\delta}{\delta-t_0}}\left[e^{-\delta\alpha^2 t}\,\mathrm{erfcx}\!\big(\alpha\sqrt{(\delta-t_0)t}\big) - e^{-t_0\alpha^2 t}\right], & t_0 < \delta,\\[2ex]
-\dfrac{2c}{\sqrt\pi}\sqrt{\dfrac{\delta}{t_0-\delta}}\; e^{-\delta\alpha^2 t}\,\mathcal{D}\!\big(\alpha\sqrt{(t_0-\delta)t}\big), & t_0 > \delta,
\end{cases}
```

where $`\mathrm{erfcx}(x) = e^{x^2}\mathrm{erfc}(x)`$ and $`\mathcal D`$ is the Dawson function (analytic
continuation of the erf branch). When both terms use the erfcx branch, the $`e^{-t_0\alpha^2 t}`$ parts
cancel exactly and are omitted — essential in water, where $`t_0<0`$ and $`e^{-t_0\alpha^2 t}`$ would
overflow. The Dawson branch had the wrong sign in the original `tllam.c`; the form in the code is
checked against the numerical inversion of the Laplace solution (`tests/validate.py`).

### 3.2 Shen TL (classic) — Baesso, Shen & Snook (1994)

Closed expression, first order in $`\theta`$ (`th`), with $`t_c = w^2/(4D)`$:

```math
I(t) = \left[1 - \frac{\theta}{2}\arctan\!\left(\frac{2mV}{\big[(1+2m)^2+V^2\big]\dfrac{t_c}{2t} + 1 + 2m + V^2}\right)\right]^2 .
```

It follows from the diffraction integral with the kernel $`K = -\theta\,\gamma(\alpha)(1-e^{-D\alpha^2 t})/\alpha`$
and $`e^{i\Phi}\approx 1+i\Phi`$. For $`|\theta|\gtrsim 0.1`$ rad the error exceeds ~1 % of the signal;
in that case the integral models are the right ones. Since there is no numerical integral, the fit
uses PSO + LM directly (`phase_based = False`).

### 3.3 TL with photobleaching — Herculano et al. (2013)

The absorption decays from $`\beta`$ to $`\bar\varepsilon\beta`$ at rate $`K_T`$ (`kt`, `eb`):
$`\beta(t) = \beta[(1-\bar\varepsilon)e^{-K_T t} + \bar\varepsilon]`$. Convolving this source with the
thermal response gives

```math
K = -\theta\,\gamma(\alpha)\frac{F(\alpha,t)}{\alpha},\qquad
F = \bar\varepsilon\,(1-e^{-bt}) + (1-\bar\varepsilon)\, b\,\frac{e^{-bt}-e^{-K_T t}}{K_T - b}.
```

$`K_T\to 0`$ or $`\bar\varepsilon = 1`$ recover the Shen TL (now without the first-order
approximation). The quotient is computed as $`t\,e^{-\min(K_T,b)t}\,(1-e^{-z})/z`$ with $`z = |K_T-b|t`$
(`_phi1`), stable when $`b\approx K_T`$.

### 3.4 TL + population lens — Silva et al. (2013)

$`\Phi = \Phi_\text{th} + \Phi_\text{PL}`$. The thermal part has the Shen kernel,
$`K = -\theta_\text{th}\gamma(\alpha)(1-e^{-D\alpha^2t})/\alpha`$. The population lens comes from the
excited-state population (two-level rate equation, saturable, lifetime $`\tau`$, excitation rate $`h`$
at the centre). It is **local** (it does not diffuse), so it depends on $`g`$ directly through the
excitation intensity $`e^{-2mg}`$:

```math
\Phi_\text{PL} = \theta_\text{PL}\int_0^l \frac{u}{1/\tau + h u}\left[1-e^{-t(1/\tau+hu)}\right]dz,
\qquad u = e^{-2mg}\,e^{-A_e z}.
```

With $`s = 1/\tau + hu`$, the integral in $`z`$ is closed:

```math
\Phi_\text{PL} = \frac{\theta_\text{PL}}{A_e h}\left[\ln\frac{s_0}{s_L} - \big(E_1(t s_L) - E_1(t s_0)\big)\right],
```

with $`E_1`$ the exponential integral ($`s_0`$ at $`z=0`$, $`s_L`$ at $`z=l`$). For $`t s_0\le 1`$ the series
$`\sum_j (-1)^{j+1} t^j(s_0^j - s_L^j)/(j\cdot j!)`$ is used (no cancellation) and, for $`A_e l<0.02`$, an
8-point Gauss–Legendre quadrature in $`z`$. As in the paper, the phase used is
$`\Phi_\text{PL}(g)-\Phi_\text{PL}(0)`$. This term is linear in $`\theta_\text{PL}`$ (`tpl`) but does not
have the form $`K\cdot[J_0-1]`$: it is evaluated directly on the $`g`$ mesh (`direct_phase`).

### 3.5 Thermal mirror LAM, finite LAM and BAM — Sato et al. (2008), Astrath et al. (2011)

The three share the form (`_TMBase`)

```math
K = 4\,\theta_\text{TM}\,\gamma(\alpha)\,F(\alpha,t),\qquad \theta_\text{TM} = -\frac{P_e A_e\alpha_T(1+\nu)\varphi}{\lambda_p k}\ [\text{1/m}],
```

and differ only in $`F`$:

| Model | $`F(\alpha,t)`$ |
|---|---|
| TM LAM (semi-infinite) | $`(1-e^{-D\alpha^2t})/\alpha^2`$ |
| TM LAM, finite sample | $`h(l\alpha)\,(1-e^{-D\alpha^2t})/\alpha^2`$ |
| TM BAM | $`\alpha D\,\hat G(\alpha,t)`$, $`\hat G = \int_0^\infty \hat T(\alpha,z,t)e^{-\alpha z}dz`$ |

In BAM, $`\hat G`$ comes from the equation

```math
\frac{d\hat G}{dt} = \frac{A_e\,\mathrm{erfc}(\alpha s) - \alpha\, e^{-\alpha^2s^2}\,\mathrm{erfcx}(A_e s)}{A_e^2-\alpha^2},\qquad s=\sqrt{Dt},
```

integrated in closed form (`_bam_closed`):

```math
\hat G = \frac{A_e I_1 - \alpha I_2}{A_e^2-\alpha^2},\quad
I_1 = \frac{2t}{x^2}\int_0^{x} x'\,\mathrm{erfc}(x')dx',\quad
I_2 = \frac{e^{-x^2}\mathrm{erfcx}(y)-1}{(A_e^2-\alpha^2)D} + \frac{A_e\,\mathrm{erf}(x)}{(A_e^2-\alpha^2)D\alpha},
```

with $`x = \alpha s`$ and $`y = A_e s`$. Where this form loses precision another one is used (see 4.6).
Limits: $`A_e\to0`$ gives LAM and $`A_e\to\infty`$ gives HAM with $`\theta_\text{HAM} = \theta_\text{TM}/A_e`$.

### 3.6 TM+TL-Fluid, HAM on/off (approximation) — Sato et al. (2008); Flizikowski et al. (2020)

Heat deposited at the surface (opaque sample), free surface. With the laser on
($`t\le t_\text{off}`$), and $`\theta = \theta_\text{HAM}`$ dimensionless:

```math
K_\text{on} = 2\theta\,\gamma(\alpha)\left[-2\sqrt{\frac{Dt}{\pi}}\,e^{-Dt\alpha^2} + \frac{\mathrm{erf}(\alpha\sqrt{Dt})}{\alpha} + 2Dt\alpha\,\mathrm{erfc}(\alpha\sqrt{Dt})\right],
```

and, after switch-off, the closed difference between $`t`$ and $`t-t_\text{off}`$. The lens in air
(active if `dndt` ≠ 0) receives heat from the surface but does not give it back:

```math
K_\text{air} = -\theta\,\frac{2\,(dn/dT)}{\beta(1+\nu)}\sqrt{\frac{D_f}{\pi}}\;\gamma(\alpha)\;
2\sqrt t\int_{\varphi_0}^{\pi/2}\cos\varphi\; e^{-\alpha^2 D_f t\sin^2\varphi}\;\mathrm{erf}\!\big(\alpha\sqrt{Dt}\cos\varphi\big)\,d\varphi ,
```

which is the integral in $`\tau`$ of the original program after the change $`\tau = t\sin^2\varphi`$
($`\varphi_0 = 0`$ on; $`\varphi_0 = \arcsin\sqrt{(t-t_\text{off})/t}`$ off). The prefactor divides
$`\theta`$ by $`\beta(1+\nu)`$, so that the lens in air is also proportional to $`\theta`$ (a single
linear parameter).

### 3.7 TM+TL-Fluid, HAM (full coupling) — Lukasievicz et al. (2012); Flizikowski et al. (2020)

Surface source, sample in $`z>0`$ and fluid in $`z<0`$, with continuity of temperature and flux at
the interface. In Hankel–Laplace space ($`q_i^2 = \alpha^2 + s/D_i`$) the solution is algebraic:

```math
\hat T(\alpha,0,s) \propto \frac{1}{s\,(k_s q_s + k_f q_f)},\quad
\hat G_s \propto \frac{1}{s\,(k_s q_s + k_f q_f)(q_s+\alpha)}\ \text{(mirror)},\quad
\hat G_f \propto \frac{1}{s\,q_f\,(k_s q_s + k_f q_f)}\ \text{(lens in the fluid)}.
```

With $`\hat s = s/(D_s\alpha^2)`$, $`\tau = D_s\alpha^2 t`$, $`\kappa = k_f/k_s`$ and $`r = D_s/D_f`$, both
responses depend only on $`\tau`$ (and on $`\kappa`$, $`r`$):

```math
g_s(\tau) = \mathcal L^{-1}\!\left\{\frac{1}{\hat s\,\big(\sqrt{1+\hat s} + \kappa\sqrt{1+r\hat s}\big)\big(\sqrt{1+\hat s}+1\big)}\right\},\qquad
g_f(\tau) = \mathcal L^{-1}\!\left\{\frac{1}{\hat s\,\sqrt{1+r\hat s}\,\big(\sqrt{1+\hat s} + \kappa\sqrt{1+r\hat s}\big)}\right\},
```

```math
K = \frac{\gamma(\alpha)}{\alpha}\left[4\,h_\text{TM}\,g_s(\tau) + \frac{h_f}{k_s}\,g_f(\tau)\right].
```

With $`\kappa=0`$ (no fluid) $`K`$ reduces to that of the TM HAM. The difference from the approximation
in 3.6 is the heat flowing back from the fluid to the sample, which enters through the term
$`\kappa\sqrt{1+r\hat s}`$.

---

## 4. Techniques that speed up the integrations

The naive computation is an integral in $`\alpha`$ (the phase) inside an integral in $`g`$ (the
diffraction), for each time, and in some models an integral in $`\tau`$ inside the one in
$`\alpha`$ as well. In the original C code each level was an adaptive quadrature (QAG) and the
integral in $`g`$ was done twice (cos and sin). The techniques below, all in `tlfit/models.py`,
replace that.

### 4.1 Separation of variables: the phase becomes a matrix product

The Bessel factor does not depend on $`t`$ or on the thermal parameters, and $`K`$ does not depend on
$`g`$. With fixed quadratures in $`\alpha`$ (nodes $`\alpha_j`$, weights $`w_j`$) and in $`g`$ (nodes
$`g_k`$):

```math
\Phi[t_i, g_k] = \sum_j K(\alpha_j, t_i)\; \underbrace{w_j\,[J_0(\alpha_j w\sqrt{m g_k}) - 1]}_{M[j,k]}
\quad\Longrightarrow\quad \Phi = K\cdot M .
```

$`M`$ depends only on the geometry ($`m`$, $`w_{0e}`$, $`V`$, accuracy, subdivision of $`g`$) and is
cached (`_bessel_matrix`, `lru_cache`): in a fit it is computed once, and each model evaluation
costs the kernel $`K`$ (in closed form) plus a matrix product (BLAS).

### 4.2 A single complex sum for the diffraction

The amplitude is $`A(t_i) = \sum_k \tilde w_k\, e^{i\varphi_{ik}}`$, with $`\varphi = Vg + \Phi`$ and the
factor $`e^{-g}`$ built into the weights ($`\tilde w_k = w_k e^{-g_k}`$). The real and imaginary parts
come out of the same pass (`cos(phi) @ wg`, `sin(phi) @ wg`), and $`I = (1+V^2)(A_r^2 + A_i^2)`$.

### 4.3 Mesh in $`\alpha`$: panelled Gauss–Legendre, in the variable $`u = \alpha w_{0e}`$

- Cut-off at $`u_\text{max} = 14`$, where $`e^{-u^2/8}\approx 2\times10^{-11}`$ (`Accuracy.u_max`).
- Panels one period of $`J_0(u\sqrt{m g_\text{max}})`$ wide (the worst case), 12 Gauss–Legendre
  nodes per panel.
- 14 levels of geometric panels near $`u=0`$ (widths $`2^{-14},\dots,2^{-1}`$ of the first panel),
  where the kernels vary as $`1/\alpha`$ and $`(1-e^{-D\alpha^2t})/\alpha^2`$.

### 4.4 Adaptive, reused mesh in $`g`$

- Base panels one period of $`e^{iVg}`$ wide (or 1, if $`V<2\pi`$), up to $`g_\text{max}=30`$.
- Each base panel is split into $`2, 4, \dots, 64`$ sub-panels where the total phase $`Vg + \Phi`$
  turns more than $`\sim2\pi`$ per panel (`needed_split`), which only happens with large $`\theta`$.
  The powers of 2 make the meshes repeat and the cache of $`M`$ get reused.
- The subdivision used in the last evaluation with the same geometry is stored (`_SPLIT_CACHE`):
  the next evaluation starts from it and almost never needs refining (`solve_phase`).

Self-convergence (comparison with doubled resolution) is ~1e-9 or better (`tests/validate.py`).

### 4.5 Integrals in $`\tau`$ solved analytically

In almost all models the time integral (convolution of the source with the thermal response) is
done in closed form in the kernel: TLLAM (erf, erfcx, Dawson), photobleaching, TM LAM, finite, BAM
and HAM, and the integral in $`z`$ of the population lens (via $`E_1`$). So the kernel is an
expression vectorized in $`(t,\alpha)`$, with no inner quadrature.

The exception is the lens in air of TM+TL-Fluid on/off (3.6), where the integral in $`\tau`$ is
numerical, but:

- the change $`\tau = t\sin^2\varphi`$ removes the $`\tau^{-1/2}`$ singularity at the end point;
- the quadrature is tanh-sinh (`tanh_sinh`, 80 nodes), good for boundary layers at the end points;
- the integral stays **outside** the loop in $`g`$ (in the original it was redone for each $`g`$,
  although it does not depend on $`g`$) and is vectorized in (time × $`\alpha`$ × $`\varphi`$), in
  blocks of 16 times to limit memory (`TMHAMTL._air`).

### 4.6 Numerically stable forms

| Where | Problem | Solution |
|---|---|---|
| $`h(x) = (\cosh x-1)/(x+\sinh x)`$ | overflow for large $`x`$ | form with `expm1(-x)` (`_hh`) |
| $`1 - e^{-x}`$ | cancellation for small $`x`$ | `-expm1(-x)` |
| TLLAM, fluid with $`t_0<0`$ | $`e^{-t_0\alpha^2t}`$ overflows | `erfcx` and analytic cancellation of the $`e^{-t_0\alpha^2 t}`$ terms; Dawson branch for $`t_0>\delta`$ |
| Photobleaching, $`b\approx K_T`$ | $`0/0`$ | $`t\,e^{-\min(K_T,b)t}(1-e^{-z})/z`$ (`_phi1`) |
| TM BAM, $`s\cdot\max(\alpha,A_e)\le 0.75`$ | cancellation in the closed form | series in $`s`$ in which $`(A_e^2-\alpha^2)`$ cancels term by term; coefficients computed once per distinct $`\alpha`$ (`_bam_series`) |
| TM BAM, $`\alpha\approx A_e`$ | removable singularity | linear interpolation between $`A_e(1\pm10^{-4})`$ |
| TM BAM, $`\int_0^X x\,\mathrm{erfc}\,x\,dx`$ with $`X<0.3`$ | cancellation | Taylor series (`_xerfc_int`) |
| Population lens, $`t s_0\le1`$ | cancellation between $`E_1`$ | alternating series |

### 4.7 Inverse Laplace tables (TM+TL-Fluid, full coupling)

$`g_s(\tau)`$ and $`g_f(\tau)`$ have no closed form. Instead of inverting the transform for each
$`(\alpha,t)`$, the code:

- inverts numerically by **fixed-contour Talbot** (Abate–Valkó, $`M=22`$, ~1e-12; `_talbot`) at 2600
  points of $`\log_{10}\tau\in[-12,\ 9.5]`$;
- interpolates by cubic spline in $`\log_{10}\tau`$ and uses the asymptotes outside the table
  ($`g_s\approx\tau/(1+\kappa\sqrt r)`$ and $`g_f\approx\tau/[\sqrt r(1+\kappa\sqrt r)]`$ for small
  $`\tau`$; $`1/[2(1+\kappa)]`$ and $`1/(1+\kappa)`$ for large $`\tau`$);
- caches the tables per $`(\kappa, r)`$ (up to ~48), so each kernel evaluation is only an
  interpolation. Since $`\kappa`$ and $`r`$ depend on $`D_s`$ and $`D_f`$, a new table is built when
  those parameters change.

### 4.8 Switch-off by superposition

The off regime reuses the same function: $`K(t) - K(t-t_\text{off})`$ (TM LAM, finite, BAM and full
coupling). In HAM on/off the difference is written in closed form.

---

## 5. Techniques in the fit

Implemented in [`tlfit/optim.py`](../tlfit/optim.py) and in `PhaseCache` (`models.py`); the
description for the user is in the README, under [Fitting method](../README.md#fitting-method).

- **Phase linear in the $`\theta`$.** For fixed non-linear parameters (typically $`D`$),
  $`\Phi = \Phi_0 + \sum_j\theta_j\Phi_j`$, with $`\Phi_j = K_j\cdot M`$ and $`K_j`$ the kernel with
  $`\theta_j = 1`$. The $`\Phi_j`$ are computed once (`PhaseCache`); varying the $`\theta`$ and I0 costs
  only the cos/sin and the sum in $`g`$, 10–15× less than a full evaluation.
- **Analytic Jacobian** in the $`\theta`$: with $`A = \sum_k\tilde w_k e^{i\varphi_k}`$,
  $`\partial A/\partial\theta_j = i\sum_k \tilde w_k\Phi_{j,k}e^{i\varphi_k}`$ and
  $`\partial I/\partial\theta_j = 2(1+V^2)\,[A_r\,\partial_j A_r + A_i\,\partial_j A_i]`$; $`\partial I/\partial I_0 = I`$.
- **Refining $`g`$ without recomputing the phase:** the maximum turning rate per panel is bounded
  by $`V + R_0 + \sum_j|\theta_j|R_j`$, with $`R_j`$ stored for each $`\Phi_j`$; the mesh is rebuilt only
  when that bound requires it.
- **Two-level fit:** sweep of the non-linear parameters (40 log-spaced points for one, 12×12 for
  two), inner fit of the $`\theta`$ at each point, and a final Levenberg–Marquardt on the best local
  minima. With more than two non-linear parameters, PSO + LM in runs.
- **Lighter search accuracy** (`Accuracy.search`): $`g_\text{max} = 18`$, meshes 0.7× and fewer
  nodes in $`\tau`$; error in $`I(t)\lesssim3\times10^{-7}`$, below the typical noise (~1e-4). The
  final refinement and the uncertainties always use normal accuracy.
- **Parallelism** with `ThreadPoolExecutor` in the sweep and in the PSO (numpy releases the GIL in
  matrix products and special functions).
