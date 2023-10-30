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
#include "types.h"
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

int solve_sub_U(RBR_INT n, RBR_INT p, double *a, RBR_INT *ia, RBR_INT *iPos) {
    std::iota(ia, ia + n, 0);

    /* find the position of nth_element */  
    std::nth_element(ia, ia + p - 1, ia + n,
            [&a](int i1, int i2) { return *(a + i1) < *(a + i2); }
            );

    /* discard those positive or zero entries */
    double *a_tmp = (double*)std::calloc(p, sizeof(double));
    RBR_INT *ia_tmp = (RBR_INT*)std::calloc(p, sizeof(RBR_INT));

    if (a_tmp == nullptr || ia_tmp == nullptr){
        std::fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return 1;
    }

    double min_pos = 0;
    RBR_INT imin_pos = -1, num = 0;
    for (RBR_INT i = 0; i < p; ++i){
        if (a[ia[i]] < 0){
            a_tmp[num] = a[ia[i]];
            ia_tmp[num] = ia[i];
            ++num;
        }
        else if (num == 0 && (imin_pos == -1 || a[ia[i]] < min_pos)){
            min_pos = a[ia[i]];
            imin_pos = ia[i];
        }
    }

    /* a can be safely discarded */
    if (num == 0){ num = 1; a[0] = -1; ia[0] = imin_pos; }
    else {
        std::memcpy(a, a_tmp, num * sizeof(double));
        std::memcpy(ia, ia_tmp, num * sizeof(int));
    }

    std::free(a_tmp);
    std::free(ia_tmp);

    *iPos = num;
    return 0;
}

int solve_sub_maxcut(RBR_INT n, RBR_INT p, double *a, RBR_INT *ia){
    //for (int i = 0; i < n; ++i) index[i] = i;
    std::iota(ia, ia + n, 0);

    /* find the position of nth_element */  
    std::nth_element(ia, ia + p - 1, ia + n,
            [&a](int i1, int i2) { return std::abs(*(a + i1)) > std::abs(*(a + i2)); }
            );

    /* rearrange */
    double *a_tmp = (double*)std::malloc(p * sizeof(double));

    if (a_tmp == nullptr){
        std::fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return 1;
    }

    for (RBR_INT i = 0; i < p; ++i){
        a_tmp[i] = a[ia[i]];
    }
    std::memcpy(a, a_tmp, p * sizeof(double));
    std::free (a_tmp);

    return 0;
}


