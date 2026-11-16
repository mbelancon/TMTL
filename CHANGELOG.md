# Notas de versão

Rascunho das notas de cada versão, para copiar na página de Releases ao criar a tag.

## Próxima versão (ainda sem número)

- **Nomenclatura:** "TM+TM-Fluid" corrigido para **TM+TL-Fluid** na interface, no README e nos
  testes. O fluido contribui sempre com uma lente térmica; não há espelho térmico no fluido (#7, PR #8).
- Ordem dos modelos no código igual à da interface e da tabela do README: TL primeiro, depois TM (#7).
- README: novo título (TMTL FIT); explica que L, B e H são as iniciais de *Low*, *Beer* e *High*
  *absorption model*; remove a seção de erros corrigidos do código original (#7).
- Documentação: novo `docs/METODOS_MATEMATICOS.md` (problema geral, fase de cada modelo, técnicas
  numéricas) e `CONTEXTO.md` na raiz do repositório (PR #9).

## 2.2.1 — errata

As notas publicadas da v2.2.1 dizem "TM+TM-Fluid, HAM on/off documentado como aproximação". O nome
correto é **TM+TL-Fluid, HAM on/off**. Texto corrigido, para substituir na página da release:

> Nomes dos modelos padronizados (TL, TM, +TL-Fluid), com a lista agrupada em TL e TM.
> TM+TL-Fluid, HAM on/off documentado como aproximação; o acoplamento completo segue Flizikowski et al. (2020) e tem teste de validação.
> Referências e agradecimentos no README.
> Corrige o nome do arquivo: a v2.2.0 saiu com "2.1.0" no nome e no LEIA-ME interno.
