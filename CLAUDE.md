# CLAUDE.md — TMTL2 (pacote tlfit)

Regras estáveis para trabalhar neste repositório. O estado atual do projeto (estrutura, decisões,
pendências) fica em `CONTEXTO.md`, importado abaixo; o README é para usuários.

@CONTEXTO.md

## Regras

- **Nomes dos modelos:** TL = lente térmica, TM = espelho térmico. Sólido + fluido (ar, água):
  **TL+TL-Fluid** ou **TM+TL-Fluid**. O fluido só contribui com TL; nunca escrever "TM-Fluid".
  LAM/BAM/HAM = Low/Beer/High Absorption Model (Sato et al., 2008).
- **Modelo novo** só entra com:
  1. referência com a equação implementada (tabela de Modelos e Referências dos dois READMEs);
  2. exemplo em `tlfit/examples/` com a linha `# tlfit: model=… m=… w0=…`;
  3. seção nos dois docs de métodos (`docs/MATH_METHODS.md` e `docs/METODOS_MATEMATICOS.md`);
  4. traduções dos textos em `tlfit/i18n.py`;
  5. teste em `tests/validate.py` (limite conhecido, forma fechada ou ajuste de dados sintéticos).
- **Idiomas:** código, comentários, mensagens e docs em inglês. Toda mudança em README ou docs vai
  para as duas línguas (`README.md` + `README.pt-BR.md`, `MATH_METHODS.md` + `METODOS_MATEMATICOS.md`).
  A interface é bilíngue: texto novo no HTML entra no dicionário `I18N` (chaves `en` e `pt`) de
  `tlfit/static/index.html`; mensagem nova no Python é `Msg(...)` com tradução em `i18n.PT`.
- **Equações no GitHub:** inline `` $`...`$ `` e blocos ` ```math `; `$...$` simples quebra com `\,`, `_` etc.
- **Numérica:** mudança em `models.py` ou `optim.py` exige `python tests/validate.py` terminando em OK,
  sem piorar os desvios contra `tests/ref/`.
- **Versão portátil:** as versões de Python, numpy e scipy são fixas em `ferramentas/montar_portatil.ps1`;
  o CI (`.github/workflows/build.yml`) monta e testa no Windows a cada push.
- **Notas de versão:** registrar cada mudança visível em `CHANGELOG.md`, na seção da próxima versão.

## Ao terminar uma sessão

Se algo mudou, propor ao Marcos a atualização do `CONTEXTO.md` como um diff curto e, aprovada,
copiar o arquivo para `CONTEXTO_TMTL.md` nos arquivos do projeto (os dois ficam iguais).
