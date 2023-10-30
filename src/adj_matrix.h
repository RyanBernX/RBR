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

#ifndef RBR_ADJ_MATRIX_H
#define RBR_ADJ_MATRIX_H

#include "types.h"
#include "api_macro.h"

typedef struct{
    RBR_INT *indx;
    RBR_INT *pntr;
    double *val;
    double *d;
    RBR_INT nnz;
    RBR_INT n;
    char label[9];
    char desc[73];
} adj_matrix;

#ifdef __cplusplus
extern "C" {
#endif

RBR_API adj_matrix * adj_matrix_create(RBR_INT n, RBR_INT nnz);
RBR_API adj_matrix * read_adj_matrix_csr(const char *filename, int is_one_based, int is_weighted);
RBR_API void adj_matrix_destroy(adj_matrix *mat);

#ifdef __cplusplus
}
#endif

#endif

