/* Substituto direto de tllam.c — mesmos argumentos, mesma saída.
 *   ./tllam_fast dados ds ths thf m w0 v l df cs cf dens denf [ns nf E coefexp dndts nu qpen qpar]
 * Os 8 últimos são opcionais (padrão = constantes do tllam.c original). */
#include "tlcore.h"
#include "tlcli_common.h"

int main(int argc, char *argv[])
{
    double p[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1.63, 1, 72, 7.5, 5.3, 0.23, 0.9e-3, 0.09e-3};
    double *t, *y, *I;
    int i, n, r;
    if (argc != 14 && argc != 22) {
        fprintf(stderr, "uso: %s dados ds ths thf m w0 v l df cs cf dens denf "
                        "[ns nf E coefexp dndts nu qpen qpar]\n", argv[0]);
        return 2;
    }
    for (i = 2; i < argc; i++) p[i - 2] = atof(argv[i]);
    n = read_data(argv[1], &t, &y);
    if (n <= 0) return 1;
    I = malloc(n * sizeof *I);
    if (tl_intensity(TL_MODEL_TLLAM, p, t, n, I, 1.0)) { fprintf(stderr, "parâmetros fora do domínio do modelo ou erro de memória\n"); return 1; }
    r = finish(argv[1], argv[2], argv[3], t, y, I, n);
    free(t); free(y); free(I);
    return r;
}
