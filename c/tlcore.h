/*
 * tlcore — núcleo numérico rápido para lente térmica / espelho térmico
 * resolvidos no tempo (substitui tllam.c e tmhamtl-onoff.c).
 *
 * Sem dependência de GSL: usa erf/erfc/j0 da libm (C99/POSIX), quadraturas de
 * Gauss-Legendre geradas internamente e uma implementação própria da função
 * de Dawson. Thread-safe (sem variáveis globais mutáveis).
 */
#ifndef TLCORE_H
#define TLCORE_H

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#define TL_API __declspec(dllexport)
#else
#define TL_API
#endif

enum { TL_MODEL_TLLAM = 0, TL_MODEL_TMHAMTL = 1 };

/* Ordem dos parâmetros
 *  TLLAM (20): d ths thf m w0 v l df cs cf dens denf ns nf E coefexp dndts nu qpen qpar
 *  TMHAMTL (10): d th m w0 v beta nu dndt df of
 * accuracy: 1.0 = padrão; 2.0 dobra a resolução das quadraturas.
 * Retorna 0 em caso de sucesso; 1 se os parâmetros levam o modelo para fora
 * do seu domínio (out recebe NaN); negativo em erro de memória/modelo. */
TL_API int tl_intensity(int model, const double *params, const double *t, int nt,
                        double *out, double accuracy);

TL_API double tl_dawson(double x);   /* exposta para testes */

#ifdef __cplusplus
}
#endif
#endif
