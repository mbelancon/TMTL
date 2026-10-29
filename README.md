# TMTL 2 — ajuste rápido de transientes de lente térmica

Ajuste de transientes de lente térmica e espelho térmico resolvidos no tempo, com interface
gráfica, ajuste de séries de transientes (várias potências)
e uma forma de integrar o modelo várias centenas de vezes mais rápida que a versão original.

> **Versão portátil para Windows, sem precisar instalar Python:** baixe o zip
> `TMTL2-…-portatil-windows.zip` na página
> [Releases](https://github.com/mbelancon/TMTL/releases), descompacte e dê dois cliques em
> `Abrir TMTL2.bat`. Detalhes em [Para usar](#para-usar).

> **Para quem usava o código original (`tllam.c`):** o termo do fluido do TLLAM tinha um
> erro de sinal que afeta amostras em água, polímeros e materiais de difusividade maior que
> a do ar. Ajustes de vidros e cristais comuns **em ar não são afetados**. Veja
> [Erros corrigidos no código original](#erros-corrigidos-no-código-original). O código
> original continua disponível em [`legado/`](legado/).

| Avaliação de um transiente | Original (C com GSL, 2 núcleos) | TMTL2 (Python/NumPy) |
|---|---|---|
| TLLAM, 86 pontos | 4,9 s | 0,02–0,05 s |
| TM-HAM on/off, 123 pontos, dn/dT ar = 0 | 81 s* | 0,2 s |
| TM-HAM on/off com lente no ar (dn/dT ≠ 0) | > 20 min (60 pontos, não terminou) | 0,5 s |

\* depois de corrigir o vazamento de memória; sem a correção o processo é morto por falta de memória.

## Modelos

| Modelo | Uso | Referência |
|---|---|---|
| **TLLAM** — TL amostra + fluido | lente térmica com acoplamento térmico amostra–fluido e espessura qualquer (modelo unificado) | Malacarne et al., *Appl. Spectrosc.* **65**, 99 (2011) |
| **TM-HAM on/off + TL no ar** | espelho térmico em amostra opaca, com desligamento e lente no ar | Sato et al., *J. Appl. Phys.* **104**, 053520 (2008) |
| **TM-LAM** | espelho térmico, baixa absorção (fonte uniforme em z), amostra semi-infinita | Sato et al. (2008), eq. 32 |
| **TM-LAM, amostra finita** | espelho térmico de baixa absorção em amostra de espessura *l* (faces livres) | Astrath et al., *J. Opt. Soc. Am. B* **28**, 1735 (2011), eq. 10 |
| **TM-BAM** (lei de Beer) | espelho térmico para qualquer A_e; tende ao LAM e ao HAM | Sato et al. (2008), eq. 33–34 |
| **TL de Shen** | expressão fechada clássica da lente térmica no modo descasado | Baesso, Shen & Snook, *J. Appl. Phys.* **75**, 3732 (1994), eq. 10 |
| **TL com fotobranqueamento** | lente térmica de líquido cuja absorção decai por reação fotoquímica (taxa K_T) até ε̄·β; ajusta θ, K_T e ε̄ | Herculano et al., *J. Phys. Chem. B* **117**, 1932 (2013), eqs. 3–6 |

**Cada modelo tem o seu θ.** A interface mostra a definição e a unidade do θ do modelo
escolhido logo abaixo da descrição, e o relatório e a série repetem essa definição:

| Modelo | θ ajustado | Unidade |
|---|---|---|
| TLLAM | θ amostra = ϑ·10⁻⁶, ϑ = P_e A_e φ/(k_s λ_p) (a fase é ϑ·l·χ(α,l)) | 10⁶ K/m |
| TLLAM | θ fluido = 2 P_e A_e φ (∂n_f/∂T)/λ_p | W/(K·m²) |
| TM-HAM | θ_HAM = θ_TM/A_e = −P_e α_T (1+ν) φ/(λ_p k) | adimensional |
| TM-LAM, TM finito, TM-BAM | θ_TM = −P_e A_e α_T (1+ν) φ/(λ_p k) | 1/m |
| Shen | θ = −P_e A_e L φ (ds/dT)/(k λ_p) | rad |
| TL com fotobranqueamento | θ = −P_e β L φ (dn/dT)/(k λ_p), β = absorção inicial (o patamar final vale ε̄·θ) | rad |

Os modelos de espelho térmico aceitam desligamento da excitação em `t_off` (0 = sem desligar). A convenção de sinal da fase é a mesma
do `tmhamtl-onoff.c` (fase somada a V·g na integral de difração); em relação à eq. 31 de
Sato et al. isso equivale a trocar o sinal de θ_TM, sem nenhum outro efeito no ajuste.

No modelo de fotobranqueamento, fixe D com uma medida do solvente puro (K_T = 0, ε̄ = 1) e ajuste
θ, K_T e ε̄ nos transientes da solução. Na série (vários P_e), K_T = K_m + K_pb·P_e: a reta de K_T
contra a potência dá K_m (difusão de massa) e K_pb (fotobranqueamento).

O modelo de Shen é de primeira ordem em θ: para |θ| ≳ 0,1 rad a diferença para a integral
de difração completa passa de ~1 % do sinal. Para sinais grandes, use os modelos integrais.

Os exemplos (`tlfit/examples/`) trazem uma linha `# tlfit: modelo=… m=… w0=…` que, ao abrir o
exemplo na interface, escolhe o modelo e preenche a geometria.

## Para usar

### Versão portátil para Windows: não precisa de Python

Na página **Releases** do repositório, baixe `TMTL2-<versão>-portatil-windows.zip`,
descompacte onde quiser e dê dois cliques em **`Abrir TMTL2.bat`**. A pasta já traz o
Python oficial (versão "embutível" do python.org), o numpy e o scipy: **não é preciso
instalar nada** nem ter direitos de administrador, e nada é alterado no computador.

A interface abre no navegador. Para encerrar, feche a aba (o programa sai sozinho cerca de
30 s depois) ou a janela preta. `Verificar instalacao.bat` confere os cálculos.

### A partir do código-fonte: para quem já tem Python

```bash
pip install numpy scipy          # o Anaconda/Miniconda já tem
python app.py                    # ou dois cliques em rodar.bat (acha o Anaconda sozinho)
```

Ou instale como pacote, direto do GitHub, e abra com `tlfit`:

```bash
pip install git+https://github.com/mbelancon/TMTL
tlfit
```

## Série de transientes (várias potências)

No seletor **Série**, no alto da barra lateral:

1. **Abrir pasta** (ou arrastar a pasta, ou **Escolher arquivos**). São lidos os arquivos
   `.dat`, `.txt`, `.csv` e `.tsv`, usando as duas primeiras colunas (tempo e sinal).
2. A potência de cada arquivo é lida do nome (`vidro_150mW.dat`, `amostra 0,2W.txt`,
   `P150.dat`) ou do cabeçalho (`# Potência: 150 mW`, `Power = 0.15 W`) e convertida
   para mW. A tabela mostra o valor e permite corrigir ou excluir arquivos. O nome da
   variável pode ser trocado (concentração, temperatura, tempo de oxidação etc.).
3. Se os sinais não estiverem normalizados (em volts, por exemplo), a normalização por
   I(t=0) é ligada sozinha; vale deixar **I0 (escala)** livre nesse caso.
4. **Ajustar série** ajusta todos os transientes com o mesmo modelo, os mesmos
   parâmetros fixos e o mesmo método (no modo automático, ~2 s por transiente no TLLAM).

Resultados:

* transientes e ajustes sobrepostos, coloridos pela potência (clique numa linha da tabela
  para destacar um transiente);
* qualquer parâmetro ajustado × potência, com barras de erro (1σ):
  * para os θ (proporcionais à potência): reta pela origem, com a inclinação θ/P_e, e
    reta com intercepto, que avisa quando o intercepto difere de zero; a opção
    "θ / potência" mostra θ/P de cada transiente e a média;
  * para D, I0 e os demais: média, desvio-padrão, erro da média e a inclinação com a
    potência, com aviso quando ela é significativa (aquecimento, efeitos não lineares);
* **Exportar resultados (.csv)**: uma linha por transiente com valores, incertezas, θ/P,
  χ² e R²; **Relatório da série (.txt)**: configuração, resultados e estatística.

As pastas `tlfit/examples/serie_vidro_ar` (normalizada) e `serie_vidro_volts` (em volts)
têm séries sintéticas para teste.

## Gerar a versão portátil

* **Automático, no GitHub:** a automação `.github/workflows/build.yml` monta e testa a versão
  portátil num Windows do GitHub a cada push (o zip fica em *Actions → execução →
  Artifacts*). Ao criar uma tag, por exemplo `git tag v2.0.0 && git push --tags`, o zip é
  publicado em **Releases**.
* **No seu computador (Windows, com internet):** dois cliques em `montar_portatil.bat`. O
  script `ferramentas/montar_portatil.ps1` baixa o Python embutível do python.org, instala
  numpy e scipy dentro dele (usando o pip de um Python que você já tenha só para baixar os
  pacotes), copia o programa, roda a validação numérica com o Python embutido e gera o zip
  em `%USERPROFILE%\TMTL2-portatil` (fora de pastas sincronizadas: são ~150 MB). O registro
  fica em `montar_portatil.log`.

As versões do Python, do numpy e do scipy são fixas no script, para que a versão portátil
seja reprodutível.

## Servidor compartilhado (opcional)

Também dá para hospedar a interface para várias pessoas: `python -m tlfit --host 0.0.0.0
--no-browser` numa máquina da rede, ou o `Dockerfile` (porta 7860) num serviço como o
Hugging Face Spaces. `TLFIT_MAX_FITS` (padrão 2) limita ajustes simultâneos; os demais
esperam na fila. Não há login.

## Método de ajuste

**Automático (padrão): ajuste em dois níveis.** Os θ e a escala I0 entram linearmente na
fase, Φ = θ_amostra·Φ_s + θ_fluido·Φ_f (no TM, Φ = θ·Φ₁). Por isso:

1. os demais parâmetros livres (tipicamente só D; no máximo dois) são varridos numa grade:
   40 pontos em escala log para um parâmetro, 12 × 12 para dois;
2. em cada ponto, Φ_s e Φ_f são calculadas uma única vez, e os θ (e I0) são ajustados
   por dentro, partindo de vários pontos, a um custo 10–15× menor que uma avaliação
   completa do modelo, com derivadas analíticas;
3. os melhores mínimos locais da grade são refinados por Levenberg–Marquardt com todos os
   parâmetros livres, na precisão normal, e daí saem as incertezas.

Como varre toda a faixa dos parâmetros não lineares, o resultado não depende de sorte, e
o registro informa quantos mínimos locais a varredura encontrou. O painel "Busca" mostra os
pontos da varredura (por exemplo, θ ótimo × D, com a cor dando o χ²). Com mais de dois
parâmetros livres além dos θ e de I0, o modo automático usa o enxame:

* **PSO** (enxame de partículas, inércia 0,9 → 0,4, topologia global) para a busca global,
  com a estimativa inicial do usuário como uma das partículas e as demais em hipercubo latino.
  Parâmetros positivos cujo intervalo cobre mais de duas décadas (difusividades) são explorados
  em escala logarítmica. O par PSO + LM é repetido em **rodadas** (padrão: 3) com sementes
  diferentes, e fica o melhor; o registro avisa se as rodadas discordam.
* **Levenberg–Marquardt** (região de confiança com limites, `scipy.optimize.least_squares`)
  a partir do melhor ponto. Converge em poucas avaliações e fornece o jacobiano,
  de onde saem as incertezas (1σ, escaladas pelo χ² reduzido) e a matriz de correlação.
  O registro avisa quando um parâmetro termina no limite ou quando dois parâmetros têm
  correlação acima de 0,95 (por exemplo, θ amostra × θ fluido no exemplo em ar).

A busca (varredura ou enxame) usa uma malha numérica mais leve: g até 18 em vez de 30 e
painéis 0,7×, com erro em I(t) abaixo de 3e-7, contra um ruído típico de ~1e-4. O refinamento
final usa sempre a precisão normal.

Tempos de um ajuste completo (2 núcleos):

| Exemplo | PSO + LM, 3 rodadas | Automático |
|---|---|---|
| TLLAM, amostra em ar | 4,7 s | 2,0 s |
| TLLAM, amostra em água | 4,9 s | 1,8 s |
| TM on/off | 22 s | 5,9 s |
| TM on/off + lente no ar, I0 livre | ~50 s (estimado) | 16 s |

Os dois modos chegam aos mesmos parâmetros nesses exemplos.

O que o notebook original fazia era uma busca aleatória em torno do melhor ponto
(perturbação de ±percent %), não um PSO: não havia velocidade nem memória das partículas,
e o alcance da busca ficava limitado pelo palpite inicial.

## O que mudou no cálculo

O modelo é

    I(t) = (1+V²) | ∫₀^30 e^{-g} exp(i[V g + Φ(g,t)]) dg |²
    Φ(g,t) = ∫ K(α,t) [J₀(α w₀ √(m g)) − 1] dα

e no código original cada Φ era uma integral adaptativa (QAG) chamada **dentro** de outra
integral adaptativa, duas vezes (uma para cos, outra para sin), e com a lente no ar havia um
terceiro nível (integral em τ) dentro da integral em α — também refeito para cada g, embora
não dependa de g.

As mudanças:

1. **Separação de variáveis.** O fator de Bessel não depende de t nem dos parâmetros
   térmicos. Com quadraturas fixas em α e g, Φ vira um produto de matrizes
   `Φ = K(t,α) · M(α,g)`, e `M` é calculada uma vez por geometria (m, w₀).
2. **Uma soma complexa** dá a parte real e a imaginária juntas.
3. **Integral em τ fora do laço em g**, com a substituição τ = t sen²φ (remove a singularidade
   τ^{-1/2}) e quadratura tanh-sinh.
4. **Quadraturas Gauss–Legendre por painéis**, com painéis geométricos perto de α = 0 e
   refinamento automático em g onde a fase total gira rápido (θ grande).
   A autoconvergência (resolução 2×) fica em ~1e-9 ou melhor.
5. Formas estáveis: `(cosh x − 1)/(x + sinh x)` sem overflow, `expm1` em `1 − e^{−x}`.

### Erros corrigidos no código original

* **TLLAM, termo do fluido — sinal trocado no ramo com a função de Dawson** (`argfasear`,
  ramos `temp0 > d` e `temp0 > df`). Com e = k/√D a efusividade térmica, a constante da
  fórmula é `temp0 = (k² − k_f²)/(e_s² − e_f²)`, e o ramo de Dawson é usado quando
  **(e_s − e_f)(D − D_f) > 0**, isto é, quando o mesmo meio tem a maior efusividade e a maior
  difusividade. Exemplos afetados: vidro em água com D acima de ~4,6e-7 m²/s, polímeros
  (PMMA) em água, amostras de difusividade maior que a do ar (silício, metais) em ar.
  **Vidros e cristais comuns em ar não são afetados** (D < D_ar). A correção sai da
  continuação analítica do ramo erf (erf(iy)/i = erfi(y)) e foi conferida contra a inversão
  numérica (Talbot) da solução exata no espaço de Laplace, em 200 casos aleatórios nos dois
  regimes (erro < 1e-7). No exemplo da amostra em água, o χ² do melhor ajuste cai de 4,3e-5
  (fórmula original) para 9,1e-6.
* TLLAM, termo do fluido — **overflow** de `exp(−temp0·α²t)` quando `temp0 < 0` (amostra
  com efusividade pouco menor que a do fluido). Os termos que crescem se cancelam
  exatamente entre si; a versão nova usa essa identidade com `erfcx` e não estoura.
* `chi` era usado sem inicialização (`reduction(+:chi)` partindo de lixo de memória).
* `tmhamtl-onoff.c`: com `dndt == 0` o workspace da GSL em `fasear` era alocado e nunca
  liberado — uma vez por avaliação do integrando — e o processo esgotava a memória.
* Workspaces da GSL eram alocados/liberados dentro dos integrandos (milhões de vezes).
* Com o manipulador de erro padrão da GSL, uma falha de arredondamento do QAG abortava o
  programa no meio do ajuste.

## Validação

`python tests/validate.py` compara com curvas geradas pelos programas originais (em
`tests/ref/`) e checa a autoconvergência:

    tllam    |novo−original| 4.9e-07   autoconvergência 2.3e-09
    tllam    |novo−original| 5.0e-07   autoconvergência 2.8e-09
    tmhamtl  |novo−original| 1.5e-04   autoconvergência 4.3e-14

Esses casos de referência são em ar, onde o código original está correto. O teste também
compara o termo do fluido com a inversão numérica da solução de Laplace nos dois regimes
da fórmula (erf e Dawson):

    fluido   Laplace, ramo erf: erro relativo máx 2.3e-08
    fluido   Laplace, ramo Dawson: erro relativo máx 3.6e-11
    ajuste   água, modo automático: χ² = 9.05310e-06, D = 4.97742e-07, θs = 4.6501, θf = -1197.1

Para os modelos de espelho térmico novos e o de Shen, o teste confere:

    TM-BAM   Ĝ(α,t) contra quadratura em τ: erro relativo máx 3.4e-08
    TM-LAM   fase contra a eq. 32 de Sato et al. (2008): erro relativo máx 2.4e-15
    limites  BAM→HAM 1.7e-09   BAM→LAM 2.1e-08   finito→LAM 1.3e-10 (relativo ao sinal)
    TL Shen  contra a integral de difração (θ = 0,01): diferença 1.1e-03 do sinal (ordem θ)
    TL branq fase contra a eq. 4 de Herculano et al. (2013): erro relativo máx 8.1e-12
    limites  K_T→0 8.8e-04   ε̄=1 8.8e-04 (relativo ao sinal; Shen é de ordem θ)
    ajuste   tmlam     D = 6.0163e-07 (+1.1σ)  θ = -120.02 (-0.3σ)
    ajuste   tmbam     D = 6.0266e-07 (+2.0σ)  θ = -149.9 (+2.1σ)
    ajuste   tmfinito  D = 5.0145e-07 (+1.4σ)  θ = -149.93 (+1.0σ)
    ajuste   tlshen    D = 4.9652e-07 (-1.3σ)  θ = -0.10419 (-1.3σ)
    ajuste   tlbranq   th = 0.048085 (+1.8σ)  kt = 2.0064 (+1.1σ)  eb = 0.27033 (+0.8σ)

* O núcleo do BAM, Ĝ(α,t) = ∫T̂(α,z,t)e^{−αz}dz, sai de
  dĜ/dt = [A_e erfc(αs) − α e^{−α²s²} erfcx(A_e s)]/(A_e² − α²), s = √(Dt), integrada em
  forma fechada. Onde a forma fechada perde precisão (s·max(α, A_e) pequeno, ou α ≈ A_e) é
  usada uma série em s em que o fator (A_e² − α²) cancela termo a termo. O resultado é
  comparado com a quadratura numérica direta de dĜ/dt em 400 pontos aleatórios
  (0,1 < A_e < 10⁸ m⁻¹, 10⁻⁶ < t < 10 s, inclusive α ≈ A_e).
* O LAM coincide com a forma fechada em g da eq. 32 de Sato et al.; o BAM tende ao HAM do
  `tmhamtl-onoff.c` (inclusive no desligamento) quando A_e cresce, e ao LAM quando A_e → 0.
* Ajustando com o LAM curvas geradas pelo BAM (geometria da fig. 11 de Sato et al.), os
  erros em θ_TM e t_c são de ~1,5 % para A_e = 77 m⁻¹ e de 16–18 % para 1061 m⁻¹, como
  no artigo: acima de ~200 m⁻¹ use o BAM.
* Os ajustes dos quatro exemplos sintéticos recuperam D e θ dentro de ±2,1σ.

No TLLAM a diferença é o arredondamento da saída do original (6 casas). No TM a diferença
vem da tolerância do original (QAG com erro relativo 1e-4): recompilando o programa original
com tolerância 1e-8, os pontos de maior discrepância (t = 0,37 s e 1,85 s) passam a coincidir
com o código novo em todas as 6 casas impressas (0,790799 e 0,998527).

## Uso como biblioteca

```python
import numpy as np
from tlfit.models import MODELS
from tlfit.optim import FitProblem, make_free, run_fit

t, y = np.loadtxt("datatofit.dat", unpack=True)
m = MODELS["tllam"]
p = m.defaults(); p.update(d=5e-7, ths=4, thf=0)
prob = FitProblem(m, t, y, p, [make_free("d", 1e-8, 1e-5), make_free("ths", 0, 20)])
res = run_fit(prob)            # modo automático (dois níveis)
print(res["params"]["d"], res["sigma"]["d"])
```

## Estrutura

    app.py              ponto de entrada (abre a interface)
    rodar.bat           roda a partir do código (acha o Anaconda/Miniconda)
    montar_portatil.bat monta a versão portátil (ferramentas/montar_portatil.ps1)
    .github/workflows/  build automático da versão portátil e Releases
    legado/             código original (C com GSL), notebook e README da versão 1
    tlfit/models.py     modelos (TLLAM, TM-HAM, TM-LAM, TM finito, TM-BAM, Shen, fotobranqueamento)
    tlfit/optim.py      ajuste em dois níveis, PSO, Levenberg–Marquardt, incertezas
    tlfit/batch.py      ajuste de séries e estatística (θ/P, médias, tendências)
    tlfit/server.py     servidor local (biblioteca padrão) para a interface
    tlfit/static/       interface (HTML/JS, sem dependências externas)
    tlfit/examples/     transientes e séries de exemplo
    tests/              validação contra o código original

Referência do modelo TLLAM: L. C. Malacarne et al., *Appl. Spectrosc.* **65**, 99 (2011),
doi:10.1366/10-06096.

## Licença

Copyright © Marcos P. Belançon e colaboradores.

Este programa é software livre, distribuído sob a **GNU General Public License, versão 3 ou
posterior** (GPL-3.0-or-later): pode ser usado, estudado, modificado e redistribuído, desde que
versões modificadas também sejam distribuídas sob a GPL. O texto completo está em
[`LICENSE`](LICENSE).

A versão portátil inclui o Python (licença PSF), o numpy e o scipy (licenças BSD); os textos
dessas licenças acompanham cada componente dentro da pasta `python` da versão portátil.
