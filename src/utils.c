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

#include <string.h>
void sparse_to_full(int n, int k, int p, const double *U, const int *iU, double *out){
    memset(out, 0, n * k * sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < p; ++j){
            int col = iU[i * p + j];
            if (col == -1) break;
            out[i * k + col] = U[i * p + j];
        }
    }
}

void sparse_to_full_c(int n, int k, int p, const double *U, const int *iU, double *out){
    memset(out, 0, n * k * sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < p; ++j){
            int col = iU[i * p + j];
            if (col == -1) break;
            out[col * n + i] = U[i * p + j];
        }
    }
}

