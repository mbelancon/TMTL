/*
 * tlcore.c — ver tlcore.h.
 *
 * Estratégia (a mesma do motor Python em tlfit/models.py):
 *
 *   I(t) = (1+V²) | ∫_0^30 e^{-g} exp(i[V g + Φ(g,t)]) dg |²
 *   Φ(g,t) = ∫ K(α,t) [J0(α w0 √(m g)) − 1] dα
 *
 * 1. O fator de Bessel não depende de t nem dos parâmetros térmicos, então a
 *    matriz M[α][g] = w_α (J0 − 1) é calculada uma vez por chamada.
 * 2. Para cada t, o núcleo K(α,t) é avaliado nos nós fixos de α e
 *    Φ(·,t) = K · M (produto vetor-matriz, em blocos de vários t).
 * 3. cos e sin saem da mesma fase: a amplitude complexa é uma única soma.
 *
 * As fórmulas dos integrandos são as do código original, sem alteração.
 */
#define _USE_MATH_DEFINES
#include "tlcore.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_MSC_VER)
#define j0 _j0
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NPP 12            /* nós de Gauss-Legendre por painel */
#define U_MAX 14.0
#define G_MAX 30.0
#define U_GRADING 14
#define MAX_SPLIT 64      /* teto de subpainéis por painel-base em g */
#define TBLOCK 8          /* tempos processados juntos no produto K·M */

/* ------------------------------------------------------------------ */
/* Funções especiais                                                   */
/* ------------------------------------------------------------------ */

/* Dawson F(x) = e^{-x²} ∫_0^x e^{s²} ds — método de Rybicki com h = 0.2
 * (erro teórico ~e^{-(π/2h)²} ≈ 1e-27) e série de Taylor perto de 0. */
double tl_dawson(double x)
{
    const double H = 0.2;
    enum { NMAX = 18 };
    static double c[NMAX];
    static int init = 0;
    double xx = fabs(x), x2, xp, e1, e2, d1, d2, sum;
    int i, n0;
    if (!init) {  /* corrida benigna: todos escrevem os mesmos valores */
        for (i = 0; i < NMAX; i++) c[i] = exp(-pow((2.0 * i + 1.0) * H, 2));
        init = 1;
    }
    if (xx < 0.05) {
        x2 = x * x;
        return x * (1 - 2.0 / 3.0 * x2 * (1 - 0.4 * x2 * (1 - 2.0 / 7.0 * x2 * (1 - 2.0 / 9.0 * x2))));
    }
    n0 = 2 * (int)floor(0.5 * xx / H + 0.5);
    xp = xx - n0 * H;
    e1 = exp(2.0 * xp * H);
    e2 = e1 * e1;
    d1 = n0 + 1;
    d2 = d1 - 2.0;
    sum = 0.0;
    for (i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
        sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
    return (x < 0 ? -1.0 : 1.0) * 0.56418958354775628695 * exp(-xp * xp) * sum;
}

/* erfcx(x) = e^{x²} erfc(x), x ≥ 0 (a libm não tem). Abaixo de 26 o produto
 * direto é exato a ~1e-13; acima, série assintótica (termos ≤ 1e-20). */
static double erfcx_pos(double x)
{
    if (x < 26.0) return exp(x * x) * erfc(x);
    {
        double x2 = x * x, s = 1.0, term = 1.0;
        int n;
        for (n = 1; n <= 6; n++) { term *= -(2.0 * n - 1.0) / (2.0 * x2); s += term; }
        return s / (x * 1.7724538509055160273);
    }
}

/* (cosh x − 1)/(x + sinh x) sem overflow */
static double hh_stable(double x)
{
    double em = expm1(-x);
    return em * em / (-expm1(-2 * x) + 2 * x * exp(-x));
}

/* ------------------------------------------------------------------ */
/* Quadraturas                                                         */
/* ------------------------------------------------------------------ */

static void gauss_legendre(int n, double *x, double *w)
{
    int i, j, m = (n + 1) / 2;
    for (i = 0; i < m; i++) {
        double z = cos(M_PI * (i + 0.75) / (n + 0.5)), z1, pp = 0;
        do {
            double p1 = 1.0, p2 = 0.0, p3;
            for (j = 0; j < n; j++) {
                p3 = p2; p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > 1e-15);
        x[i] = -z; x[n - 1 - i] = z;
        w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
}

/* nós/pesos compostos para os painéis definidos por edges[0..ne-1] */
static int composite(const double *edges, int ne, double **xo, double **wo)
{
    double gx[NPP], gw[NPP];
    int i, k, n = (ne - 1) * NPP;
    double *x = malloc(n * sizeof *x), *w = malloc(n * sizeof *w);
    if (!x || !w) { free(x); free(w); return -1; }
    gauss_legendre(NPP, gx, gw);
    for (i = 0; i < ne - 1; i++) {
        double a = edges[i], b = edges[i + 1];
        for (k = 0; k < NPP; k++) {
            x[i * NPP + k] = 0.5 * (b - a) * gx[k] + 0.5 * (b + a);
            w[i * NPP + k] = 0.5 * (b - a) * gw[k];
        }
    }
    *xo = x; *wo = w;
    return n;
}

static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static int u_grid(double m, double w0, double f, double **u, double **wu)
{
    double u_top = fmin(U_MAX, 2e5 * w0);
    double period = 2 * M_PI / sqrt(fmax(m, 1e-12) * G_MAX);
    double width = period / f, first = fmin(width, u_top);
    int nun = (int)ceil((u_top - first) / width);
    int i, ne = 0, nmax, n;
    double *e;
    if (nun < 1) nun = 1;
    nmax = 1 + U_GRADING + nun + 1;
    e = malloc(nmax * sizeof *e);
    if (!e) return -1;
    e[ne++] = 0.0;
    for (i = U_GRADING; i >= 1; i--) e[ne++] = first * pow(2.0, -i);
    for (i = 0; i <= nun; i++) e[ne++] = first + (u_top - first) * i / nun;
    qsort(e, ne, sizeof *e, cmp_dbl);
    for (i = 1, n = 1; i < ne; i++) if (e[i] != e[n - 1]) e[n++] = e[i];
    i = composite(e, n, u, wu);
    free(e);
    return i;
}

/* Painéis-base em g: largura = um período de e^{i v g} (ou 1 se v < 2π). */
static int g_base(double v, double f, double **edges)
{
    double rate = fmax(fabs(v), 2 * M_PI);
    double width = 2 * M_PI / rate / f;
    int i, n = (int)ceil(G_MAX / width);
    double *e = malloc((n + 1) * sizeof *e);
    if (!e) return -1;
    for (i = 0; i <= n; i++) e[i] = G_MAX * i / n;
    *edges = e;
    return n;  /* número de painéis */
}

/* nós em g com o painel-base i dividido em split[i] subpainéis */
static int g_grid(const double *base, int nb, const int *split, double **g, double **wg)
{
    int i, k, ne = 1, r;
    double *e;
    for (i = 0; i < nb; i++) ne += split[i];
    e = malloc(ne * sizeof *e);
    if (!e) return -1;
    e[0] = base[0];
    for (i = 0, ne = 1; i < nb; i++)
        for (k = 1; k <= split[i]; k++)
            e[ne++] = base[i] + (base[i + 1] - base[i]) * k / split[i];
    r = composite(e, ne, g, wg);
    free(e);
    return r;
}

static int tanh_sinh(int n, double *x, double *w)
{
    double h = 6.0 / n;
    int k, m = 0;
    for (k = -n / 2; k <= n / 2; k++) {
        double kk = k * h, s = 0.5 * M_PI * sinh(kk);
        double xx = 1.0 / (1.0 + exp(-2 * s));
        double ww = 0.5 * h * 0.5 * M_PI * cosh(kk) / pow(cosh(s), 2);
        if (xx > 0 && xx < 1 && ww > 1e-300) { x[m] = xx; w[m] = ww; m++; }
    }
    return m;
}

/* ------------------------------------------------------------------ */
/* Núcleos K(α, t)                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    int model;
    const double *p;
    int nu;
    const double *alpha;
    double *B;               /* TLLAM: fator de material da amostra, por α */
    int ntau;
    double tx[400], tw[400]; /* tanh-sinh para o ar (TMHAMTL) */
} Ctx;

static void kernel_tllam(const Ctx *c, double t, double *K)
{
    const double *p = c->p;
    double d = p[0], ths = p[1], thf = p[2], w0 = p[4], l = p[6], df = p[7];
    double k = d * p[8] * p[10], kf = df * p[9] * p[11];
    double k2 = k * k, kf2 = kf * kf;
    double t0 = (k2 - kf2) * d * df / (k2 * df - kf2 * d);
    int j;
    for (j = 0; j < c->nu; j++) {
        double a = c->alpha[j], a2 = a * a;
        double gauss = exp(-a2 * w0 * w0 / 8), val = 0.0;
        if (ths != 0)
            val += ths * l * gauss * c->B[j] * (-expm1(-d * t * a2)) / a;
        if (thf != 0) {
            double temp1 = gauss / ((k2 - kf2) * a2);
            double bb1 = k * erf(a * sqrt(df * t)) - kf * erf(a * sqrt(d * t));
            /* Cada termo é A·e^{−t0α²t} + R; no ramo erf do original,
             *   −c s e^{−t0α²t} erf(x) = −c s e^{−t0α²t} + c s e^{−ddα²t} erfcx(x).
             * Com os dois termos nesse ramo, k s(df) = kf s(d) e as partes
             * e^{−t0α²t} se cancelam exatamente (essencial para t0 < 0). */
            double A2 = 0, A3 = 0, R2, R3, bb;
            if (!isfinite(t0)) { /* t0 → ±∞: esses termos tendem a zero */
                val += thf * temp1 * bb1;
                K[j] = val;
                continue;
            }
            if (t0 > d)
                /* sinal corrigido em relação ao tllam.c (ver tlfit/models.py) */
                R3 = -(2 * kf / sqrt(M_PI)) * sqrt(d / (t0 - d)) * exp(-a2 * t * d) *
                     tl_dawson(a * sqrt((t0 - d) * t));
            else {
                double sc = kf * sqrt(d / (d - t0));
                A3 = -sc;
                R3 = sc * exp(-d * a2 * t) * erfcx_pos(a * sqrt((d - t0) * t));
            }
            if (t0 > df)
                R2 = -(2 * k / sqrt(M_PI)) * sqrt(df / (t0 - df)) * exp(-a2 * t * df) *
                     tl_dawson(a * sqrt((t0 - df) * t));
            else {
                double sc = k * sqrt(df / (df - t0));
                A2 = -sc;
                R2 = sc * exp(-df * a2 * t) * erfcx_pos(a * sqrt((df - t0) * t));
            }
            bb = R2 - R3;
            if (!(t0 < d && t0 < df)) bb += (A2 - A3) * exp(-t0 * a2 * t);   /* aqui t0 > 0 */
            val += thf * temp1 * (bb1 + bb);
        }
        K[j] = val;
    }
}

static void kernel_tmhamtl(const Ctx *c, double t, double *K)
{
    const double *p = c->p;
    double d = p[0], th = p[1], w0 = p[3], beta = p[5], nu = p[6], dndt = p[7];
    double df = p[8], of = p[9];
    double st = sqrt(d * t);
    double pref = 0, lo = 0;
    int j, q;
    if (dndt != 0) {
        pref = -th * 2 * dndt / (beta * (1 + nu)) * sqrt(df / M_PI);
        lo = (t <= of) ? 0.0 : asin(sqrt((t - of) / t));
    }
    for (j = 0; j < c->nu; j++) {
        double a = c->alpha[j], a2 = a * a;
        double gauss = exp(-a2 * w0 * w0 / 8), val;
        if (t <= of) {
            val = 2 * th * gauss * (-2 * sqrt(d * t / M_PI) * exp(-d * t * a2)
                                    + erf(st * a) / a + erfc(st * a) * 2 * d * t * a);
        } else {
            double to = t - of, sto = sqrt(d * to);
            val = 2 * th * (2 * sqrt(d * to / M_PI) * exp(-a2 * (d * to + w0 * w0 / 8))
                            - 2 * sqrt(d * t / M_PI) * exp(-a2 * (d * t + w0 * w0 / 8))
                            + gauss / a * ((2 * d * t * a2 - 1) * erfc(st * a)
                                           - (2 * d * to * a2 - 1) * erfc(sto * a)));
        }
        if (dndt != 0) {
            /* ∫ τ^{-1/2} e^{-α² df τ} erf(α√(d(t−τ))) dτ com τ = t sin²φ */
            double span = 0.5 * M_PI - lo, J = 0;
            for (q = 0; q < c->ntau; q++) {
                double ph = lo + span * c->tx[q], s = sin(ph), co = cos(ph);
                J += span * c->tw[q] * co * exp(-a2 * df * t * s * s) * erf(a * st * co);
            }
            val += pref * gauss * 2 * sqrt(t) * J;
        }
        K[j] = val;
    }
}

/* ------------------------------------------------------------------ */
/* Intensidade                                                         */
/* ------------------------------------------------------------------ */

static void eval_kernel(const Ctx *c, double t, double *K)
{
    if (t <= 1e-9) { memset(K, 0, c->nu * sizeof *K); return; }
    if (c->model == TL_MODEL_TLLAM) kernel_tllam(c, t, K);
    else kernel_tmhamtl(c, t, K);
}

int tl_intensity(int model, const double *params, const double *t, int nt,
                 double *out, double accuracy)
{
    double m, w0, v, *u = NULL, *wu = NULL, *g = NULL, *wg = NULL, *M = NULL;
    double *alpha = NULL, *B = NULL, *K = NULL, *phi = NULL, *base = NULL, *rmax = NULL;
    int *split = NULL, *nsplit = NULL;
    int nu, ng = 0, nb, i, j, k, err = 0, it;
    Ctx ctx;
    if (accuracy <= 0) accuracy = 1.0;
    if (model == TL_MODEL_TLLAM) { m = params[3]; w0 = params[4]; v = params[5]; }
    else if (model == TL_MODEL_TMHAMTL) { m = params[2]; w0 = params[3]; v = params[4]; }
    else return -2;

    nu = u_grid(m, w0, accuracy, &u, &wu);
    nb = g_base(v, accuracy, &base);
    if (nu < 0 || nb < 0) { err = -1; goto done; }
    alpha = malloc(nu * sizeof *alpha);
    B = malloc(nu * sizeof *B);
    K = malloc((size_t)nt * nu * sizeof *K);
    split = malloc(nb * sizeof *split);
    nsplit = malloc(nb * sizeof *nsplit);
    rmax = malloc(nb * sizeof *rmax);
    if (!alpha || !B || !K || !split || !nsplit || !rmax) { err = -1; goto done; }
    for (j = 0; j < nu; j++) alpha[j] = u[j] / w0;
    for (i = 0; i < nb; i++) split[i] = 1;

    memset(&ctx, 0, sizeof ctx);
    ctx.model = model; ctx.p = params; ctx.nu = nu; ctx.alpha = alpha; ctx.B = B;
    if (model == TL_MODEL_TLLAM) {
        const double *p = params;
        double l = p[6], ns = p[12], nf = p[13], E = p[14], ce = p[15], dndts = p[16];
        double nup = p[17], qpen = p[18], qpar = p[19];
        double qci2 = pow(ns, 3) * E * ce / (4 * (1 - nup));
        for (j = 0; j < nu; j++) {
            double x = l * alpha[j], hl = hh_stable(x) / x;
            B[j] = dndts + 4 * (ns - nf) * (1 + nup) * ce * hl
                   + qci2 * (qpar + 3 * qpen - 4 * (qpar * nup + qpen * (2 + nup)) * hl);
        }
    } else {
        int nn = (int)(80 * accuracy);
        if (nn > 390) nn = 390;
        ctx.ntau = tanh_sinh(nn, ctx.tx, ctx.tw);
    }

    /* núcleo K(α, t) para todos os tempos */
#pragma omp parallel for schedule(dynamic)
    for (i = 0; i < nt; i++) eval_kernel(&ctx, t[i], K + (size_t)i * nu);
    for (i = 0; i < nt * nu; i++)
        if (!isfinite(K[i])) {      /* parâmetros fora do domínio do modelo */
            for (i = 0; i < nt; i++) out[i] = NAN;
            err = 1;
            goto done;
        }

    for (it = 0; it < 6; it++) {
        int changed = 0;
        free(g); free(wg); free(M); free(phi);
        g = wg = M = phi = NULL;
        ng = g_grid(base, nb, split, &g, &wg);
        if (ng < 0) { err = -1; goto done; }
        M = malloc((size_t)nu * ng * sizeof *M);
        phi = malloc((size_t)nt * ng * sizeof *phi);
        if (!M || !phi) { err = -1; goto done; }
#pragma omp parallel for private(k) schedule(static)
        for (j = 0; j < nu; j++)
            for (k = 0; k < ng; k++)
                M[(size_t)j * ng + k] = wu[j] / w0 * (j0(u[j] * sqrt(m * g[k])) - 1.0);

        /* Φ = K·M em blocos de TBLOCK tempos */
#pragma omp parallel for schedule(dynamic)
        for (i = 0; i < nt; i += TBLOCK) {
            int nbk = (nt - i < TBLOCK) ? nt - i : TBLOCK, b, jj, kk;
            for (b = 0; b < nbk; b++)
                for (kk = 0; kk < ng; kk++) phi[(size_t)(i + b) * ng + kk] = v * g[kk];
            for (jj = 0; jj < nu; jj++) {
                const double *Mr = M + (size_t)jj * ng;
                for (b = 0; b < nbk; b++) {
                    double kv = K[(size_t)(i + b) * nu + jj];
                    double *ph = phi + (size_t)(i + b) * ng;
                    for (kk = 0; kk < ng; kk++) ph[kk] += kv * Mr[kk];
                }
            }
        }

        /* refinamento: fase total deve girar ≤ ~2π por painel */
        for (i = 0; i < nb; i++) rmax[i] = 0;
        for (k = 0; k + 1 < ng; k++) {
            double r = 0, mid = 0.5 * (g[k] + g[k + 1]);
            int bi = (int)(mid / (base[1] - base[0]));
            if (bi >= nb) bi = nb - 1;
            for (i = 0; i < nt; i++) {
                double dd = fabs(phi[(size_t)i * ng + k + 1] - phi[(size_t)i * ng + k]);
                if (dd > r) r = dd;
            }
            r /= (g[k + 1] - g[k]);
            if (r > rmax[bi]) rmax[bi] = r;
        }
        for (i = 0; i < nb; i++) {
            double need = rmax[i] * (base[i + 1] - base[i]) / (2 * M_PI);
            int s = 1;
            if (!(need <= MAX_SPLIT)) need = MAX_SPLIT;   /* também pega NaN */
            if (need > 1) s = 1 << (int)ceil(log2(need));
            nsplit[i] = s > split[i] ? s : split[i];
            if (nsplit[i] != split[i]) changed = 1;
        }
        if (!changed) break;
        memcpy(split, nsplit, nb * sizeof *split);
    }

    for (k = 0; k < ng; k++) wg[k] *= exp(-g[k]);
#pragma omp parallel for private(k)
    for (i = 0; i < nt; i++) {
        double re = 0, im = 0;
        const double *ph = phi + (size_t)i * ng;
        if (t[i] <= 1e-9) { out[i] = 1.0; continue; }
        for (k = 0; k < ng; k++) { re += wg[k] * cos(ph[k]); im += wg[k] * sin(ph[k]); }
        out[i] = (re * re + im * im) * (1 + v * v);
    }

done:
    free(u); free(wu); free(g); free(wg); free(M); free(alpha); free(B);
    free(K); free(phi); free(base); free(split); free(nsplit); free(rmax);
    return err;
}
