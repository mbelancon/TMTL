# CONTEXTO — TMTL2 (pacote tlfit)

Arquivo para o Claude (e para quem for mexer no código). Curto e atualizado; o README é para usuários.
Fica na raiz do repositório; a cópia `CONTEXTO_TMTL.md` nos arquivos do projeto deve acompanhar esta.

## Objetivo
Ajustar transientes de lente térmica (TL) e espelho térmico (TM) resolvidos no tempo, com interface
gráfica e ajuste de séries (vários transientes em potências diferentes → parâmetros vs. potência,
que é a análise usual do grupo). Distribuição para colaboradores via GitHub.

## Estrutura
- `tlfit/` — `models.py`, `optim.py`, `batch.py`, `dataio.py`, `server.py` + `static/` (interface), `examples/`
- `app.py`, `rodar.bat`, `montar_portatil.bat`, `Dockerfile`, `tests/`
- `docs/METODOS_MATEMATICOS.md` — física e numérica dos modelos (ler antes de mexer em `models.py`)
- `legado/` — referência histórica; o código em C foi abandonado (TMTL2 é só Python desde a v2.0.0)
- Versão atual no `pyproject.toml` (2.2.1 em nov/2026); licença GPL-3

## Métodos matemáticos (resumo; detalhes em `docs/METODOS_MATEMATICOS.md`)
- Todos os modelos integrais: I(t) = (1+V²)|∫₀^30 e^{−g} e^{i[Vg+Φ]} dg|², Φ = ∫ K(α,t)[J₀(α w0 √(mg)) − 1] dα
- Cada modelo só implementa `kernel(t, a, p, acc)`; a difração é comum. Exceções: TL de Shen (fechado) e lente de população (`direct_phase`)
- θ entram linearmente em K → ajuste em dois níveis (`PhaseCache`); desligamento = K(t) − K(t − t_off)
- Aceleração: Φ = K·M com M em cache por geometria; Gauss–Legendre por painéis; g adaptativo com subdivisão em cache; núcleos em forma fechada; tabelas de Talbot + spline no TM+TL-Fluid completo
- Modelo novo: escrever o núcleo, a seção em `docs/METODOS_MATEMATICOS.md` e um exemplo em `tlfit/examples/`

## Convenções de nomes
- TM = espelho térmico, TL = lente térmica
- Sólido + fluido (ar/água): TL+TL-Fluid ou TM+TL-Fluid. O fluido só contribui com TL; nunca existe "TM no fluido"
- LAM/BAM/HAM = Low/Beer/High Absorption Model (Sato et al., 2008)

## Decisões
- Distribuição: zip portátil Windows (Python embutido, não precisa instalar) + código-fonte; não usar PyInstaller
- `tmhamtl` mantido como aproximação sem retorno do fluido (documentado no README e na interface), ao lado do acoplamento completo
- Ajuste conjunto TL+TM deixado de lado; modelos por EDP em espera
- Versão web no site pessoal: por enquanto só divulgação, sem rodar no navegador

## Referências-chave
- Shen 10.1063/1.356045; TL + população 10.1063/1.4801838; fotobranqueamento (Herculano) 10.1021/jp3119296
- Flizikowski et al., Opt. Express 28, 7116 (2020) — TM-HAM on/off
- Artigos de referência no Google Drive (pasta Claude)

## Pendências
- (preencher)
