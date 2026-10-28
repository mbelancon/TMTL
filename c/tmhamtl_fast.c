/* Substituto direto de tmhamtl-onoff.c — mesmos argumentos, mesma saída.
 *   ./tmhamtl_fast dados d th m w0 v beta nu dndt df of */
#include "tlcore.h"
#include "tlcli_common.h"

int main(int argc, char *argv[])
{
    double p[10], *t, *y, *I;
    int i, n, r;
    if (argc != 12) {
        fprintf(stderr, "uso: %s dados d th m w0 v beta nu dndt df of\n", argv[0]);
        return 2;
    }
    for (i = 2; i < argc; i++) p[i - 2] = atof(argv[i]);
    n = read_data(argv[1], &t, &y);
    if (n <= 0) return 1;
    I = malloc(n * sizeof *I);
    if (tl_intensity(TL_MODEL_TMHAMTL, p, t, n, I, 1.0)) { fprintf(stderr, "parâmetros fora do domínio do modelo ou erro de memória\n"); return 1; }
    r = finish(argv[1], argv[2], argv[3], t, y, I, n);
    free(t); free(y); free(I);
    return r;
}
