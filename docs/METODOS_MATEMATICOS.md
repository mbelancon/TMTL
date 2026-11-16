# Métodos matemáticos

Este documento descreve o que os modelos do `tlfit` calculam, por que a **fase** é a única coisa
que muda de um modelo para outro e quais técnicas numéricas tornam o cálculo rápido. As equações
seguem o código de [`tlfit/models.py`](../tlfit/models.py) (a notação entre crases é a do código);
as referências de cada modelo estão na tabela de Modelos do [README](../README.md#modelos).

Convenção de nomes: **TL** = lente térmica, **TM** = espelho térmico, **+TL-Fluid** = lente térmica
do fluido em contato com a amostra (o fluido só contribui com TL; não há TM no fluido).
LAM/BAM/HAM = *Low/Beer/High Absorption Model*.

Sumário

1. [O problema geral](#1-o-problema-geral)
2. [A fase é o que diferencia os modelos](#2-a-fase-é-o-que-diferencia-os-modelos)
3. [A particularidade de cada modelo](#3-a-particularidade-de-cada-modelo)
4. [Técnicas para acelerar as integrações](#4-técnicas-para-acelerar-as-integrações)
5. [Técnicas no ajuste](#5-técnicas-no-ajuste)

---

## 1. O problema geral

Todos os modelos descrevem o mesmo experimento, em três etapas:

1. **Aquecimento.** O feixe de excitação, gaussiano com raio $w_{0e}$ na amostra (`w0`), deposita
   calor com perfil radial $e^{-2r^2/w_{0e}^2}$ e um perfil em profundidade $Q(z)$ que depende da
   absorção. A temperatura $T(r,z,t)$ obedece à equação de difusão

   ```math
   \frac{\partial T}{\partial t} - D\,\nabla^2 T = Q_0\, e^{-2r^2/w_{0e}^2}\, Q(z),
   ```

   com $D$ a difusividade térmica (`d`) e as condições de contorno de cada modelo (superfície
   livre, amostra finita, troca de calor com o fluido).

2. **Mudança de fase do feixe de prova.** O aquecimento altera o caminho óptico: o índice de
   refração muda com a temperatura (lente térmica, na amostra ou no fluido) e a superfície se
   deforma por dilatação (espelho térmico). O feixe de prova, de raio $w_{1p}$ na amostra, sai com
   uma fase $\Phi(r,t)$ que é medida em relação ao centro, $\Phi(r,t) = \phi(r,t) - \phi(0,t)$.

3. **Propagação até o detector.** O feixe de prova se propaga (difração de Fresnel) até o
   detector, no campo distante, que mede só o centro do feixe. Com a variável radial
   $g = (r/w_{1p})^2$, a intensidade normalizada é

   ```math
   I(t) = (1+V^2)\,\left|\int_0^{\infty} e^{-g}\, e^{\,i\,[V g + \Phi(g,t)]}\, dg\right|^2 ,
   ```

   onde $m = (w_{1p}/w_{0e})^2$ e $V$ (`m`, `v`) são os parâmetros geométricos do arranjo
   descasado ($V \approx Z_1/Z_c$, com $Z_1$ a distância da cintura do prova à amostra e $Z_c$ o
   comprimento confocal). Sem fase, $\int_0^\infty e^{-(1-iV)g}dg = 1/(1-iV)$ e o fator $(1+V^2)$
   faz $I = 1$: o sinal é normalizado pelo valor antes do aquecimento. O código integra $g$ até
   $g_\text{max} = 30$ (o termo desprezado é da ordem de $e^{-30}$), como o programa original.

A escala `scale` (I0) multiplica $I(t)$ e serve para sinais que não estão exatamente normalizados.

## 2. A fase é o que diferencia os modelos

As etapas 1 e 3 são iguais em todos os modelos (a não ser pelo $Q(z)$ e pelas fronteiras). O que
muda é **como a temperatura vira fase**:

| Efeito | Fase (em relação ao centro) |
|---|---|
| **TL** (amostra de espessura $L$) | $\Phi = \dfrac{2\pi}{\lambda_p}\dfrac{ds}{dT}\displaystyle\int_0^L [T(r,z,t)-T(0,z,t)]\,dz$ |
| **TM** (superfície refletora) | $\Phi = \dfrac{2\pi}{\lambda_p}\, 2\,[u_z(r,0,t)-u_z(0,0,t)]$, com $u_z \propto -(1+\nu)\alpha_T \displaystyle\int_0^\infty \hat T(\alpha,z,t)\,e^{-\alpha z}dz$ |
| **TL-Fluid** | $\Phi_f = \dfrac{2\pi}{\lambda_p}\dfrac{\partial n_f}{\partial T}\displaystyle\int_\text{fluido} [T_f(r,z,t)-T_f(0,z,t)]\,dz$ (dois passes no TM) |
| Lente de população | variação do índice causada pela população do estado excitado (não térmica) |

Na TL a fase vem da temperatura **integrada na espessura**; no TM vem do **deslocamento da
superfície**, que pesa a temperatura com $e^{-\alpha z}$ (solução termoelástica quase estática).
É por isso que um mesmo perfil de temperatura dá fases diferentes em TL e TM.

### Forma comum: o núcleo $K(\alpha,t)$

Escrevendo a temperatura pela transformada de Hankel, $T(r,z,t) = \int_0^\infty \hat T(\alpha,z,t)\,
J_0(\alpha r)\,\alpha\,d\alpha$, e usando $r = w_{1p}\sqrt{g} = w_{0e}\sqrt{m g}$, todas as fases
térmicas ficam com a mesma forma:

```math
\Phi(g,t) = \int_0^\infty K(\alpha,t)\,\left[J_0\!\left(\alpha\, w_{0e}\sqrt{m g}\right) - 1\right] d\alpha .
```

O fator $[J_0 - 1]$ é a parte "geométrica" (onde o prova enxerga a amostra) e é o mesmo em todos
os modelos. **Toda a física de cada modelo está no núcleo $K(\alpha,t)$**: a fonte, as fronteiras,
o efeito óptico (TL, TM, TL-Fluid) e o desligamento. A transformada de Hankel da fonte gaussiana
dá o fator $e^{-\alpha^2 w_{0e}^2/8}$ que aparece em todos os núcleos. No código, cada modelo
implementa só o método `kernel(t, a, p, acc)`; a integral de difração (`Model.intensity`) é
compartilhada.

Duas consequências práticas:

- os parâmetros $\theta$ entram **linearmente** em $K$ e, portanto, em $\Phi$:
  $\Phi = \theta_\text{amostra}\Phi_s + \theta_\text{fluido}\Phi_f$ (no TM, $\Phi = \theta\,\Phi_1$).
  Isso é usado no ajuste em dois níveis (seção 5);
- o desligamento da excitação em $t_\text{off}$ (`of`) é a superposição
  $K_\text{off}(t) = K(t) - K(t-t_\text{off})$, porque a equação do calor é linear.

Exceções: o **TL de Shen** é uma expressão fechada para $I(t)$ (não passa pela integral numérica)
e a **lente de população** tem uma fase que não tem a forma $K\cdot[J_0-1]$, calculada diretamente
em $g$ (`direct_phase`).

## 3. A particularidade de cada modelo

Notação: $a=\alpha$, $w = w_{0e}$, $\gamma(\alpha) = e^{-\alpha^2 w^2/8}$, $b = D\alpha^2$ e
$h(x) = (\cosh x - 1)/(x + \sinh x)$ (fator de espessura finita, `_hh`).

### Resumo

| Modelo | Fonte $Q(z)$ | Fronteiras | Fase | Núcleo $K(\alpha,t)$ |
|---|---|---|---|---|
| TL+TL-Fluid (TLLAM) | uniforme (baixa absorção) | espessura $l$, troca de calor com o fluido | TL amostra + TL fluido | fechado; fluido com erf, erfcx, Dawson |
| TL de Shen | uniforme | amostra fina, sem fluido | TL | — (forma fechada de $I(t)$) |
| TL com fotobranqueamento | uniforme, decaindo no tempo | amostra fina | TL | fechado |
| TL + lente de população | Beer | amostra fina | TL + população | TL fechado; população direta em $g$ |
| TM LAM | uniforme | semi-infinita, superfície livre | TM | fechado |
| TM LAM, amostra finita | uniforme | espessura $l$, faces livres | TM | fechado, com $h(l\alpha)$ |
| TM BAM | Beer, $A_e e^{-A_e z}$ | semi-infinita | TM | fechado + série perto de $\alpha s\to 0$ |
| TM+TL-Fluid, HAM on/off (aproximação) | na superfície | semi-infinita, superfície livre | TM + TL fluido (ar) | TM fechado; fluido por quadratura em $\tau$ |
| TM+TL-Fluid, HAM (acoplamento completo) | na superfície | troca de calor com o fluido | TM + TL fluido | Laplace–Hankel invertido numericamente |

### 3.1 TL+TL-Fluid (TLLAM) — Malacarne et al. (2011)

Lente térmica da amostra (espessura $l$, `l`) mais a do fluido, com acoplamento térmico. O núcleo
é a soma de dois termos lineares em $\theta_s$ (`ths`) e $\theta_f$ (`thf`):

```math
K = \theta_s\, l\,\gamma(\alpha)\,B(\alpha)\,\frac{1-e^{-D\alpha^2 t}}{\alpha} \;+\; \theta_f\,K_f(\alpha,t).
```

O fator $B(\alpha)$ (o $\chi$ do README) junta as três contribuições ópticas da amostra: o
$dn/dT$, a deformação das duas faces e o efeito fotoelástico,

```math
B = \frac{dn}{dT} + 4(n_s-n_f)(1+\nu)\alpha_T\frac{h(l\alpha)}{l\alpha}
  + \frac{n_s^3 E\alpha_T}{4(1-\nu)}\left[q_\parallel + 3q_\perp - 4\big(q_\parallel\nu + q_\perp(2+\nu)\big)\frac{h(l\alpha)}{l\alpha}\right].
```

O termo do fluido, com $k = D c_s\rho_s$, $k_f = D_f c_f \rho_f$ e
$t_0 = (k^2-k_f^2) D D_f/(k^2 D_f - k_f^2 D)$, é

```math
K_f = \frac{\gamma(\alpha)}{(k^2-k_f^2)\alpha^2}\Big[k\,\mathrm{erf}(\alpha\sqrt{D_f t}) - k_f\,\mathrm{erf}(\alpha\sqrt{D t})
      + \mathcal{T}(D_f, k) - \mathcal{T}(D, k_f)\Big],
```

com, para cada par $(\delta, c)$,

```math
\mathcal{T}(\delta,c) =
\begin{cases}
c\sqrt{\dfrac{\delta}{\delta-t_0}}\left[e^{-\delta\alpha^2 t}\,\mathrm{erfcx}\!\big(\alpha\sqrt{(\delta-t_0)t}\big) - e^{-t_0\alpha^2 t}\right], & t_0 < \delta,\\[2ex]
-\dfrac{2c}{\sqrt\pi}\sqrt{\dfrac{\delta}{t_0-\delta}}\; e^{-\delta\alpha^2 t}\,\mathcal{D}\!\big(\alpha\sqrt{(t_0-\delta)t}\big), & t_0 > \delta,
\end{cases}
```

onde $\mathrm{erfcx}(x) = e^{x^2}\mathrm{erfc}(x)$ e $\mathcal D$ é a função de Dawson (continuação
analítica do ramo erf). Quando os dois termos usam o ramo erfcx, as partes $e^{-t_0\alpha^2 t}$ se
cancelam exatamente e são omitidas — essencial em água, onde $t_0<0$ e $e^{-t_0\alpha^2 t}$
estouraria. O ramo de Dawson tinha o sinal trocado no `tllam.c` original; a forma do código é
conferida contra a inversão numérica da solução de Laplace (`tests/validate.py`).

### 3.2 TL de Shen (clássico) — Baesso, Shen & Snook (1994)

Expressão fechada, de primeira ordem em $\theta$ (`th`), com $t_c = w^2/(4D)$:

```math
I(t) = \left[1 - \frac{\theta}{2}\arctan\!\left(\frac{2mV}{\big[(1+2m)^2+V^2\big]\dfrac{t_c}{2t} + 1 + 2m + V^2}\right)\right]^2 .
```

Ela sai da integral de difração com o núcleo $K = -\theta\,\gamma(\alpha)(1-e^{-D\alpha^2 t})/\alpha$
e $e^{i\Phi}\approx 1+i\Phi$. Para $|\theta|\gtrsim 0{,}1$ rad o erro passa de ~1 % do sinal; nesse
caso os modelos integrais são os adequados. Como não há integral numérica, o ajuste usa PSO + LM
direto (`phase_based = False`).

### 3.3 TL com fotobranqueamento — Herculano et al. (2013)

A absorção decai de $\beta$ até $\bar\varepsilon\beta$ com taxa $K_T$ (`kt`, `eb`):
$\beta(t) = \beta[(1-\bar\varepsilon)e^{-K_T t} + \bar\varepsilon]$. A convolução dessa fonte com a
resposta térmica dá

```math
K = -\theta\,\gamma(\alpha)\frac{F(\alpha,t)}{\alpha},\qquad
F = \bar\varepsilon\,(1-e^{-bt}) + (1-\bar\varepsilon)\, b\,\frac{e^{-bt}-e^{-K_T t}}{K_T - b}.
```

$K_T\to 0$ ou $\bar\varepsilon = 1$ recuperam a TL de Shen (agora sem a aproximação de primeira
ordem). O quociente é calculado como $t\,e^{-\min(K_T,b)t}\,(1-e^{-z})/z$ com $z = |K_T-b|t$
(`_phi1`), estável quando $b\approx K_T$.

### 3.4 TL + lente de população — Silva et al. (2013)

$\Phi = \Phi_\text{th} + \Phi_\text{PL}$. A parte térmica tem o núcleo de Shen,
$K = -\theta_\text{th}\gamma(\alpha)(1-e^{-D\alpha^2t})/\alpha$. A lente de população vem da
população do estado excitado (equação de taxa de dois níveis, saturável, vida $\tau$, taxa de
excitação $h$ no centro). Ela é **local** (não difunde), então depende de $g$ diretamente pela
intensidade da excitação $e^{-2mg}$:

```math
\Phi_\text{PL} = \theta_\text{PL}\int_0^l \frac{u}{1/\tau + h u}\left[1-e^{-t(1/\tau+hu)}\right]dz,
\qquad u = e^{-2mg}\,e^{-A_e z}.
```

Com $s = 1/\tau + hu$, a integral em $z$ é fechada:

```math
\Phi_\text{PL} = \frac{\theta_\text{PL}}{A_e h}\left[\ln\frac{s_0}{s_L} - \big(E_1(t s_L) - E_1(t s_0)\big)\right],
```

com $E_1$ a integral exponencial ($s_0$ em $z=0$, $s_L$ em $z=l$). Para $t s_0\le 1$ usa-se a série
$\sum_j (-1)^{j+1} t^j(s_0^j - s_L^j)/(j\cdot j!)$ (sem cancelamento) e, para $A_e l<0{,}02$, uma
quadratura de Gauss–Legendre de 8 pontos em $z$. Como no artigo, a fase usada é
$\Phi_\text{PL}(g)-\Phi_\text{PL}(0)$. Esse termo é linear em $\theta_\text{PL}$ (`tpl`) mas não tem a
forma $K\cdot[J_0-1]$: é avaliado direto na malha de $g$ (`direct_phase`).

### 3.5 Espelho térmico LAM, LAM finito e BAM — Sato et al. (2008), Astrath et al. (2011)

Os três compartilham a forma (`_TMBase`)

```math
K = 4\,\theta_\text{TM}\,\gamma(\alpha)\,F(\alpha,t),\qquad \theta_\text{TM} = -\frac{P_e A_e\alpha_T(1+\nu)\varphi}{\lambda_p k}\ [\text{1/m}],
```

e diferem só em $F$:

| Modelo | $F(\alpha,t)$ |
|---|---|
| TM LAM (semi-infinita) | $(1-e^{-D\alpha^2t})/\alpha^2$ |
| TM LAM, amostra finita | $h(l\alpha)\,(1-e^{-D\alpha^2t})/\alpha^2$ |
| TM BAM | $\alpha D\,\hat G(\alpha,t)$, $\hat G = \int_0^\infty \hat T(\alpha,z,t)e^{-\alpha z}dz$ |

No BAM, $\hat G$ sai da equação

```math
\frac{d\hat G}{dt} = \frac{A_e\,\mathrm{erfc}(\alpha s) - \alpha\, e^{-\alpha^2s^2}\,\mathrm{erfcx}(A_e s)}{A_e^2-\alpha^2},\qquad s=\sqrt{Dt},
```

integrada em forma fechada (`_bam_closed`):

```math
\hat G = \frac{A_e I_1 - \alpha I_2}{A_e^2-\alpha^2},\quad
I_1 = \frac{2t}{x^2}\int_0^{x} x'\,\mathrm{erfc}(x')dx',\quad
I_2 = \frac{e^{-x^2}\mathrm{erfcx}(y)-1}{(A_e^2-\alpha^2)D} + \frac{A_e\,\mathrm{erf}(x)}{(A_e^2-\alpha^2)D\alpha},
```

com $x = \alpha s$ e $y = A_e s$. Onde essa forma perde precisão usa-se outra (ver 4.6). Limites:
$A_e\to0$ dá o LAM e $A_e\to\infty$ dá o HAM com $\theta_\text{HAM} = \theta_\text{TM}/A_e$.

### 3.6 TM+TL-Fluid, HAM on/off (aproximação) — Sato et al. (2008); Flizikowski et al. (2020)

Calor depositado na superfície (amostra opaca), superfície livre. Com o laser ligado
($t\le t_\text{off}$), $\theta = \theta_\text{HAM}$ adimensional:

```math
K_\text{on} = 2\theta\,\gamma(\alpha)\left[-2\sqrt{\frac{Dt}{\pi}}\,e^{-Dt\alpha^2} + \frac{\mathrm{erf}(\alpha\sqrt{Dt})}{\alpha} + 2Dt\alpha\,\mathrm{erfc}(\alpha\sqrt{Dt})\right],
```

e, depois do desligamento, a diferença fechada entre $t$ e $t-t_\text{off}$. A lente no ar (ativa se
`dndt` ≠ 0) recebe o calor da superfície mas não o devolve:

```math
K_\text{ar} = -\theta\,\frac{2\,(dn/dT)}{\beta(1+\nu)}\sqrt{\frac{D_f}{\pi}}\;\gamma(\alpha)\;
2\sqrt t\int_{\varphi_0}^{\pi/2}\cos\varphi\; e^{-\alpha^2 D_f t\sin^2\varphi}\;\mathrm{erf}\!\big(\alpha\sqrt{Dt}\cos\varphi\big)\,d\varphi ,
```

que é a integral em $\tau$ do programa original depois da troca $\tau = t\sin^2\varphi$
($\varphi_0 = 0$ ligado; $\varphi_0 = \arcsin\sqrt{(t-t_\text{off})/t}$ desligado). O prefator divide
$\theta$ por $\beta(1+\nu)$, de modo que a lente do ar também é proporcional a $\theta$ (um único
parâmetro linear).

### 3.7 TM+TL-Fluid, HAM (acoplamento completo) — Lukasievicz et al. (2012); Flizikowski et al. (2020)

Fonte superficial, amostra em $z>0$ e fluido em $z<0$, com continuidade de temperatura e fluxo
na interface. No espaço de Hankel–Laplace ($q_i^2 = \alpha^2 + s/D_i$) a solução é algébrica:

```math
\hat T(\alpha,0,s) \propto \frac{1}{s\,(k_s q_s + k_f q_f)},\quad
\hat G_s \propto \frac{1}{s\,(k_s q_s + k_f q_f)(q_s+\alpha)}\ \text{(espelho)},\quad
\hat G_f \propto \frac{1}{s\,q_f\,(k_s q_s + k_f q_f)}\ \text{(lente no fluido)}.
```

Com $\hat s = s/(D_s\alpha^2)$, $\tau = D_s\alpha^2 t$, $\kappa = k_f/k_s$ e $r = D_s/D_f$, as duas
respostas dependem só de $\tau$ (e de $\kappa$, $r$):

```math
g_s(\tau) = \mathcal L^{-1}\!\left\{\frac{1}{\hat s\,\big(\sqrt{1+\hat s} + \kappa\sqrt{1+r\hat s}\big)\big(\sqrt{1+\hat s}+1\big)}\right\},\qquad
g_f(\tau) = \mathcal L^{-1}\!\left\{\frac{1}{\hat s\,\sqrt{1+r\hat s}\,\big(\sqrt{1+\hat s} + \kappa\sqrt{1+r\hat s}\big)}\right\},
```

```math
K = \frac{\gamma(\alpha)}{\alpha}\left[4\,h_\text{TM}\,g_s(\tau) + \frac{h_f}{k_s}\,g_f(\tau)\right].
```

Com $\kappa=0$ (sem fluido) $K$ se reduz ao do TM HAM. A diferença para a aproximação 3.6 é o
retorno de calor do fluido para a amostra, que entra pelo termo $\kappa\sqrt{1+r\hat s}$.

---

## 4. Técnicas para acelerar as integrações

O cálculo ingênuo é uma integral em $\alpha$ (a fase) dentro de uma integral em $g$ (a difração),
para cada tempo, e em alguns modelos ainda uma integral em $\tau$ dentro da de $\alpha$. No código
original em C cada nível era uma quadratura adaptativa (QAG) e a integral em $g$ era feita duas
vezes (cos e sin). As técnicas abaixo, todas em `tlfit/models.py`, substituem isso.

### 4.1 Separação de variáveis: a fase vira um produto de matrizes

O fator de Bessel não depende de $t$ nem dos parâmetros térmicos, e $K$ não depende de $g$. Com
quadraturas fixas em $\alpha$ (nós $\alpha_j$, pesos $w_j$) e em $g$ (nós $g_k$):

```math
\Phi[t_i, g_k] = \sum_j K(\alpha_j, t_i)\; \underbrace{w_j\,[J_0(\alpha_j w\sqrt{m g_k}) - 1]}_{M[j,k]}
\quad\Longrightarrow\quad \Phi = K\cdot M .
```

$M$ depende só da geometria ($m$, $w_{0e}$, $V$, precisão, subdivisão de $g$) e é guardada em
cache (`_bessel_matrix`, `lru_cache`): num ajuste ela é calculada uma vez e cada avaliação do
modelo custa o núcleo $K$ (em forma fechada) mais um produto de matrizes (BLAS).

### 4.2 Uma soma complexa para a difração

A amplitude é $A(t_i) = \sum_k \tilde w_k\, e^{i\varphi_{ik}}$, com $\varphi = Vg + \Phi$ e o fator
$e^{-g}$ embutido nos pesos ($\tilde w_k = w_k e^{-g_k}$). A parte real e a imaginária saem da mesma
passada (`cos(phi) @ wg`, `sin(phi) @ wg`), e $I = (1+V^2)(A_r^2 + A_i^2)$.

### 4.3 Malha em $\alpha$: Gauss–Legendre por painéis, na variável $u = \alpha w_{0e}$

- Corte em $u_\text{max} = 14$, onde $e^{-u^2/8}\approx 2\times10^{-11}$ (`Accuracy.u_max`).
- Painéis com largura de um período de $J_0(u\sqrt{m g_\text{max}})$ (o pior caso), 12 nós de
  Gauss–Legendre por painel.
- 14 níveis de painéis geométricos perto de $u=0$ (larguras $2^{-14},\dots,2^{-1}$ do primeiro
  painel), onde os núcleos variam como $1/\alpha$ e $(1-e^{-D\alpha^2t})/\alpha^2$.

### 4.4 Malha em $g$ adaptativa e reaproveitada

- Painéis-base com largura de um período de $e^{iVg}$ (ou 1, se $V<2\pi$), até $g_\text{max}=30$.
- Cada painel-base é subdividido em $2, 4, \dots, 64$ subpainéis onde a fase total
  $Vg + \Phi$ gira mais de $\sim2\pi$ por painel (`needed_split`), o que só acontece com $\theta$
  grande. As potências de 2 fazem as malhas se repetirem e o cache de $M$ ser reaproveitado.
- A subdivisão usada na última avaliação com a mesma geometria fica guardada (`_SPLIT_CACHE`):
  a próxima avaliação começa por ela e quase sempre não precisa refinar (`solve_phase`).

A autoconvergência (comparação com a resolução dobrada) fica em ~1e-9 ou melhor
(`tests/validate.py`).

### 4.5 Integrais em $\tau$ resolvidas analiticamente

Em quase todos os modelos a integral no tempo (convolução da fonte com a resposta térmica) é
feita em forma fechada no núcleo: TLLAM (erf, erfcx, Dawson), fotobranqueamento, TM LAM,
finito, BAM e HAM, e a integral em $z$ da lente de população (via $E_1$). Assim o núcleo é uma
expressão vetorizada em $(t,\alpha)$, sem quadratura interna.

A exceção é a lente no ar do TM+TL-Fluid on/off (3.6), em que a integral em $\tau$ é numérica,
mas:

- a troca $\tau = t\sin^2\varphi$ remove a singularidade $\tau^{-1/2}$ na ponta;
- a quadratura é tanh-sinh (`tanh_sinh`, 80 nós), boa para camadas-limite nas pontas;
- a integral fica **fora** do laço em $g$ (no original era refeita para cada $g$, embora não
  dependa de $g$) e é vetorizada em (tempo × $\alpha$ × $\varphi$), em blocos de 16 tempos para
  limitar a memória (`TMHAMTL._air`).

### 4.6 Formas numericamente estáveis

| Onde | Problema | Solução |
|---|---|---|
| $h(x) = (\cosh x-1)/(x+\sinh x)$ | overflow para $x$ grande | forma com `expm1(-x)` (`_hh`) |
| $1 - e^{-x}$ | cancelamento para $x$ pequeno | `-expm1(-x)` |
| TLLAM, fluido com $t_0<0$ | $e^{-t_0\alpha^2t}$ estoura | `erfcx` e cancelamento analítico dos termos $e^{-t_0\alpha^2 t}$; ramo de Dawson para $t_0>\delta$ |
| Fotobranqueamento, $b\approx K_T$ | $0/0$ | $t\,e^{-\min(K_T,b)t}(1-e^{-z})/z$ (`_phi1`) |
| TM BAM, $s\cdot\max(\alpha,A_e)\le 0{,}75$ | cancelamento na forma fechada | série em $s$ em que $(A_e^2-\alpha^2)$ cancela termo a termo; coeficientes calculados uma vez por $\alpha$ distinto (`_bam_series`) |
| TM BAM, $\alpha\approx A_e$ | singularidade removível | interpolação linear entre $A_e(1\pm10^{-4})$ |
| TM BAM, $\int_0^X x\,\mathrm{erfc}\,x\,dx$ com $X<0{,}3$ | cancelamento | série de Taylor (`_xerfc_int`) |
| Lente de população, $t s_0\le1$ | cancelamento entre $E_1$ | série alternada |

### 4.7 Tabelas da inversa de Laplace (TM+TL-Fluid, acoplamento completo)

$g_s(\tau)$ e $g_f(\tau)$ não têm forma fechada. Em vez de inverter a transformada para cada
$(\alpha,t)$, o código:

- inverte numericamente por **Talbot de contorno fixo** (Abate–Valkó, $M=22$, ~1e-12; `_talbot`)
  em 2600 pontos de $\log_{10}\tau\in[-12;\ 9{,}5]$;
- interpola por spline cúbica em $\log_{10}\tau$ e usa as assíntotas fora da tabela
  ($g_s\approx\tau/(1+\kappa\sqrt r)$ e $g_f\approx\tau/[\sqrt r(1+\kappa\sqrt r)]$ para $\tau$
  pequeno; $1/[2(1+\kappa)]$ e $1/(1+\kappa)$ para $\tau$ grande);
- guarda as tabelas em cache por $(\kappa, r)$ (até ~48), então cada avaliação do núcleo é só
  uma interpolação. Como $\kappa$ e $r$ dependem de $D_s$ e $D_f$, uma tabela nova é feita quando
  esses parâmetros mudam.

### 4.8 Desligamento por superposição

O regime desligado reaproveita a mesma função: $K(t) - K(t-t_\text{off})$ (TM LAM, finito, BAM e
acoplamento completo). No HAM on/off a diferença é escrita em forma fechada.

---

## 5. Técnicas no ajuste

Implementadas em [`tlfit/optim.py`](../tlfit/optim.py) e em `PhaseCache` (`models.py`); a descrição
para o usuário está no README, em [Método de ajuste](../README.md#método-de-ajuste).

- **Fase linear nos $\theta$.** Para os parâmetros não lineares fixos (tipicamente $D$),
  $\Phi = \Phi_0 + \sum_j\theta_j\Phi_j$, com $\Phi_j = K_j\cdot M$ e $K_j$ o núcleo com $\theta_j = 1$.
  Os $\Phi_j$ são calculados uma vez (`PhaseCache`); variar os $\theta$ e I0 custa só o cos/sin e a
  soma em $g$, 10–15× menos que uma avaliação completa.
- **Jacobiano analítico** nos $\theta$: com $A = \sum_k\tilde w_k e^{i\varphi_k}$,
  $\partial A/\partial\theta_j = i\sum_k \tilde w_k\Phi_{j,k}e^{i\varphi_k}$ e
  $\partial I/\partial\theta_j = 2(1+V^2)\,[A_r\,\partial_j A_r + A_i\,\partial_j A_i]$; $\partial I/\partial I_0 = I$.
- **Refinamento de $g$ sem recalcular a fase:** a taxa máxima de giro por painel é limitada por
  $V + R_0 + \sum_j|\theta_j|R_j$, com $R_j$ guardado para cada $\Phi_j$; a malha só é refeita
  quando esse limite exige.
- **Ajuste em dois níveis:** varredura dos não lineares (40 pontos em log para um, 12×12 para
  dois), ajuste interno dos $\theta$ em cada ponto, e Levenberg–Marquardt final nos melhores
  mínimos locais. Com mais de dois não lineares, PSO + LM em rodadas.
- **Precisão de busca mais leve** (`Accuracy.search`): $g_\text{max} = 18$, malhas 0,7× e menos nós
  em $\tau$; erro em $I(t)\lesssim3\times10^{-7}$, abaixo do ruído típico (~1e-4). O refinamento final
  e as incertezas usam sempre a precisão normal.
- **Paralelismo** por `ThreadPoolExecutor` na varredura e no PSO (o numpy libera o GIL nos
  produtos de matrizes e nas funções especiais).
