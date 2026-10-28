/* Leitura de dados e saída compartilhadas pelos programas de linha de comando. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Lê duas colunas (tempo, intensidade); ignora linhas vazias e comentários (#). */
static int read_data(const char *path, double **t, double **y)
{
    FILE *f = fopen(path, "r");
    char line[512];
    int n = 0, cap = 256;
    if (!f) { fprintf(stderr, "não foi possível abrir %s\n", path); return -1; }
    *t = malloc(cap * sizeof **t); *y = malloc(cap * sizeof **y);
    while (fgets(line, sizeof line, f)) {
        double a, b;
        char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\r' || *s == 0) continue;
        for (char *c = s; *c; c++) if (*c == ',' || *c == ';') *c = ' ';
        if (sscanf(s, "%lf %lf", &a, &b) != 2) continue;
        if (n == cap) { cap *= 2; *t = realloc(*t, cap * sizeof **t); *y = realloc(*y, cap * sizeof **y); }
        (*t)[n] = a; (*y)[n] = b; n++;
    }
    fclose(f);
    return n;
}

static int finish(const char *datafile, const char *a2, const char *a3,
                  const double *t, const double *y, const double *I, int n)
{
    char path[1024];
    double chi = 0.0;
    int i;
    FILE *f;
    for (i = 0; i < n; i++) chi += (y[i] - I[i]) * (y[i] - I[i]);
    printf("%.10lf\n", chi);
    snprintf(path, sizeof path, "result-%s-%s-%s", a2, a3, datafile);
    f = fopen(path, "w");
    if (!f) { fprintf(stderr, "não foi possível criar %s\n", path); return 1; }
    for (i = 0; i < n; i++) fprintf(f, "%.10g\t%.10g\n", t[i], I[i]);
    fclose(f);
    return 0;
}
