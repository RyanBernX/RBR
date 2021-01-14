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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_MKL
#include "mkl.h"
#else
#include "cblas.h"
#endif
#include "DCRBR.h"

void cal_confusion_matrix(int n, int k, const int *labels, const int *info, int *mat){
    memset(mat, 0, k * k * sizeof(int));

    /* loop over every node */
    for (int i = 0; i < n; ++i){
        /* C(ci, cj) */
        int ci = labels[i];
        int cj = info[i];

        ++mat[cj * k + ci];
    }
}

double cal_Fsame(int n, int k, const int *labels, const int *info){
    double Fsame = 0;
    int *C = (int*)malloc(k * k * sizeof(double));

    /* obtain confustion matrix */
    cal_confusion_matrix(n, k, labels, info, C);

    for (int i = 0; i < k; ++i){
        double max = 0;
        for (int j = 0; j < k; ++j)
            if (max < C[j * k + i]) max = C[j * k + i];

        Fsame += max;
    }

    for (int i = 0; i < k; ++i){
        double max = 0;
        for (int j = 0; j < k; ++j)
            if (max < C[i * k + j]) max = C[i * k + j];

        Fsame += max;
    }

    return Fsame / (2 * n);
}

double cal_NMI(int n, int k, const int *labels, const int *info){
    int Nt = 0, *Ni, *Nj, *C;

    Ni = (int*)calloc(k, sizeof(int));
    Nj = (int*)calloc(k, sizeof(int));
    C  = (int*)malloc(k * k * sizeof(double));

    /* obtain confustion matrix */
    cal_confusion_matrix(n, k, labels, info, C);
    /* compute Ni, Nj, Nt */
    for (int i = 0; i < k; ++i){
        for (int j = 0; j < k; ++j){
            Ni[i] += C[j * k + i];
            Nj[i] += C[i * k + j];
        }
        Nt += Ni[i];
    }
    double nom = 0, den = 0;
    for (int i = 0; i < k; ++i){
        for (int j = 0; j < k; ++j){
            double Nij = C[j * k + i];
            nom += Nij * log((Nij * Nt) / (Ni[i] * Nj[j]));
        }
    }
    for (int i = 0; i < k; ++i){
        den += Ni[i] * log(1.0 * Ni[i] / Nt) + Nj[i] * log(1.0 * Nj[i] / Nt);
    }

    return -2 * nom / den;
}

int cal_mis(int n, int k, const int *labels, const int *info) {
    int *cnt, mis = 0;
    cnt = (int*)calloc(k, sizeof(int));

    for (int i = 0; i < k; ++i) {
        int cntU = 0;
        for (int j = 0; j < n; ++j)
            if (labels[j] == i) {
                ++cntU;
                cnt[info[j]]++;
            }

        int m = cnt[0];
        for (int j = 1; j < k; ++j)
            if (m < cnt[j]) m = cnt[j];

        mis += cntU - m;
        memset(cnt, 0, sizeof(int)*k);
    }
    free(cnt);
    return mis;
}


/** @fn double cal_Q()
 * @brief compute the modularity Q
 */
double cal_Q(adj_matrix *A, int n, int k, const int *labels, double lambda, const double *UTB){
    double obj = 0;
    int *ir, *jc;

    ir = A->pntr; jc = A->indx;

    /* compute tr(U' * A * U) */
    for (int i = 0; i < n; ++i){
        for (int j = ir[i]; j < ir[i + 1]; ++j){
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
double cal_Strength(adj_matrix *A, int n, int k, const int *labels){
    double *scores, ret;
    int *deg, *pntr = A->pntr, *indx = A->indx;

    scores = (double*)malloc(k * sizeof(double));
    deg = (int*)calloc(k, sizeof(int));

    for (int i = 0; i < k; ++i) scores[i] = 1;

    for (int i = 0; i < n; ++i){ /* for each node */
        int deg_in = 0, deg_out = 0;
        for (int j = pntr[i]; j < pntr[i+1]; ++j){
            int t = indx[j]; /* destnation */
            labels[i] == labels[t] ? ++deg_in : ++deg_out;
        }
        if (deg_in <= deg_out) scores[labels[i]] = 0.5;

        deg[labels[i]] += deg_in - deg_out;
    }

    for (int i = 0; i < k; ++i)
        if (deg[i] <= 0) scores[i] = 0;

    ret = cblas_dasum(k, scores, 1) / k;
    free(scores); free(deg);
    return ret;
}

double cal_CC(adj_matrix *A, int n, int k, const int *labels){
    double *scores, ret;
    int *comm_cnt;

    scores = (double*)calloc(k, sizeof(double));
    comm_cnt = (int*)calloc(k, sizeof(int));

    for (int v = 0; v < n; ++v){ /* for each node */
        int label = labels[v], deg = A->d[v], CC_cnt = 0;
        int *s = A->indx + A->pntr[v];
        ++comm_cnt[label];

        if (deg <= 1){
            scores[label] += 1;
            continue;
        }

        int dom = deg * (deg - 1) / 2;
        /* enumerate all neighbours */
        for (int i = 0; i < deg; ++i)
            for (int  j = i + 1; j < deg; ++j)
                if (labels[s[i]] == label && labels[s[j]] == label && have_edge(A, s[i], s[j]))
                    ++CC_cnt;

        scores[label] += 1.0 * CC_cnt / dom;
    }

    for (int i = 0; i < k; ++i)
        scores[i] = comm_cnt[i] > 0 ? scores[i] / comm_cnt[i] : 1;

    ret = cblas_dasum(k, scores, 1) / k;
    free(scores); free(comm_cnt);
    return ret;
}

int have_edge(adj_matrix *A, int s, int t){
    int *indx = A->indx, *pntr = A->pntr;

    for (int j = pntr[s]; j < pntr[s+1]; ++j){
        /* TODO: use binary search if it is too slow */
        if (indx[j] == t) return 1;
    }

    return 0;
}
