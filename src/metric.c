/*
 * RBR
 *
 * Copyright (C) 2019  Haoyang Liu (liuhaoyang@pku.edu.cn)
 *                     Zaiwen Wen  (wenzw@pku.edu.cn)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.\
 */

#include "rbr_subroutines.h"
#include "rbr_blas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void cal_confusion_matrix(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info, RBR_INT *mat){
    memset(mat, 0, k * k * sizeof(RBR_INT));

    /* loop over every node */
    for (RBR_INT i = 0; i < n; ++i){
        /* C(ci, cj) */
        RBR_INT ci = labels[i];
        RBR_INT cj = info[i];

        ++mat[cj * k + ci];
    }
}

double cal_Fsame(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info){
    double Fsame = 0;
    RBR_INT *C = (RBR_INT*)malloc(k * k * sizeof(RBR_INT));
    
    if (C == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return NAN;
    }

    /* obtain confustion matrix */
    cal_confusion_matrix(n, k, labels, info, C);

    for (RBR_INT i = 0; i < k; ++i){
        double max = 0;
        for (RBR_INT j = 0; j < k; ++j)
            if (max < C[j * k + i]) max = C[j * k + i];

        Fsame += max;
    }

    for (RBR_INT i = 0; i < k; ++i){
        double max = 0;
        for (RBR_INT j = 0; j < k; ++j)
            if (max < C[i * k + j]) max = C[i * k + j];

        Fsame += max;
    }

    return Fsame / (2 * n);
}

double cal_NMI(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info){
    RBR_INT Nt = 0, *Ni, *Nj, *C;

    Ni = (RBR_INT*)calloc(k, sizeof(RBR_INT));
    Nj = (RBR_INT*)calloc(k, sizeof(RBR_INT));
    C  = (RBR_INT*)malloc(k * k * sizeof(RBR_INT));

    if (Ni == NULL || Nj == NULL || C == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        if (Ni) free(Ni);
        if (Nj) free(Nj);
        if (C) free(C);
        return NAN;
    }

    /* obtain confustion matrix */
    cal_confusion_matrix(n, k, labels, info, C);
    /* compute Ni, Nj, Nt */
    for (RBR_INT i = 0; i < k; ++i){
        for (RBR_INT j = 0; j < k; ++j){
            Ni[i] += C[j * k + i];
            Nj[i] += C[i * k + j];
        }
        Nt += Ni[i];
    }
    double nom = 0, den = 0;
    for (RBR_INT i = 0; i < k; ++i){
        for (RBR_INT j = 0; j < k; ++j){
            double Nij = C[j * k + i];
            nom += Nij * log((Nij * Nt) / (Ni[i] * Nj[j]));
        }
    }
    for (RBR_INT i = 0; i < k; ++i){
        den += Ni[i] * log(1.0 * Ni[i] / Nt) + Nj[i] * log(1.0 * Nj[i] / Nt);
    }

    return -2 * nom / den;
}

RBR_INT cal_mis(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info) {
    RBR_INT *cnt, mis = 0;
    cnt = (RBR_INT*)calloc(k, sizeof(RBR_INT));

    if (cnt == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return -1;
    }

    for (RBR_INT i = 0; i < k; ++i) {
        RBR_INT cntU = 0;
        for (RBR_INT j = 0; j < n; ++j)
            if (labels[j] == i) {
                ++cntU;
                cnt[info[j]]++;
            }

        RBR_INT m = cnt[0];
        for (RBR_INT j = 1; j < k; ++j)
            if (m < cnt[j]) m = cnt[j];

        mis += cntU - m;
        memset(cnt, 0, sizeof(RBR_INT)*k);
    }
    free(cnt);
    return mis;
}


/** @fn double cal_Q()
 * @brief compute the modularity Q
 */
double cal_Q(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels, double lambda, const double *UTB){
    double obj = 0;
    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;

    /* compute tr(U' * A * U) */
    for (RBR_INT i = 0; i < n; ++i){
        for (RBR_INT j = ir[i]; j < ir[i + 1]; ++j){
            /* at (i, jc[j]) */
            if (labels[i] == labels[jc[j]]) obj += 1;
        }
    }

    obj -= lambda * cblas_ddot(k, UTB, 1, UTB, 1);

    return obj / A->nnz;
}


/** @fn double cal_Strength()
 * @brief compute the strength based on the results
 */
double cal_Strength(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels){
    double *scores, ret;

    const RBR_INT *pntr = A->pntr;
    const RBR_INT *indx = A->indx;

    scores = (double*)malloc(k * sizeof(double));
    RBR_INT *deg = (RBR_INT*)calloc(k, sizeof(RBR_INT));

    if (scores == NULL || deg == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        if (scores) free(scores);
        if (deg) free(deg);

        return NAN;
    }

    for (RBR_INT i = 0; i < k; ++i) scores[i] = 1;

    for (RBR_INT i = 0; i < n; ++i){ /* for each node */
        RBR_INT deg_in = 0, deg_out = 0;
        for (RBR_INT j = pntr[i]; j < pntr[i+1]; ++j){
            RBR_INT t = indx[j]; /* destnation */
            labels[i] == labels[t] ? ++deg_in : ++deg_out;
        }
        if (deg_in <= deg_out) scores[labels[i]] = 0.5;

        deg[labels[i]] += deg_in - deg_out;
    }

    for (RBR_INT i = 0; i < k; ++i)
        if (deg[i] <= 0) scores[i] = 0;

    ret = cblas_dasum(k, scores, 1) / k;
    free(scores);
    free(deg);
    return ret;
}

double cal_CC(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels){
    double *scores, ret;
    RBR_INT *comm_cnt;

    scores = (double*)calloc(k, sizeof(double));
    comm_cnt = (RBR_INT*)calloc(k, sizeof(RBR_INT));

    if (scores == NULL || comm_cnt == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        if (scores) free(scores);
        if (comm_cnt) free(comm_cnt);

        return NAN;
    }

    for (RBR_INT v = 0; v < n; ++v){ /* for each node */
        RBR_INT label = labels[v], deg = A->d[v], CC_cnt = 0;
        RBR_INT *s = A->indx + A->pntr[v];
        ++comm_cnt[label];

        if (deg <= 1){
            scores[label] += 1;
            continue;
        }

        RBR_INT dom = deg * (deg - 1) / 2;
        /* enumerate all neighbours */
        for (RBR_INT i = 0; i < deg; ++i)
            for (RBR_INT  j = i + 1; j < deg; ++j)
                if (labels[s[i]] == label && labels[s[j]] == label && have_edge(A, s[i], s[j]))
                    ++CC_cnt;

        scores[label] += 1.0 * CC_cnt / dom;
    }

    for (RBR_INT i = 0; i < k; ++i)
        scores[i] = comm_cnt[i] > 0 ? scores[i] / comm_cnt[i] : 1;

    ret = cblas_dasum(k, scores, 1) / k;
    free(scores);
    free(comm_cnt);
    return ret;
}

int have_edge(const adj_matrix *A, RBR_INT s, RBR_INT t){
    const RBR_INT *indx = A->indx, *pntr = A->pntr;

    for (RBR_INT j = pntr[s]; j < pntr[s+1]; ++j){
        /* TODO: use binary search if it is too slow */
        if (indx[j] == t) return 1;
    }

    return 0;
}
