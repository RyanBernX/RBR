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

#ifndef RBR_RBR_H
#define RBR_RBR_H

#include "adj_matrix.h"
#include "types.h"
#include "api_macro.h"

#ifdef __cplusplus
extern "C" {
#endif

// main functions
RBR_API rbr_result * rbr(const adj_matrix *A, RBR_INT k, const rbr_param *param);
RBR_API rbr_result * rbr_maxcut(const adj_matrix *A, RBR_INT k, const rbr_param *param);

// rbr_results
// defined in rbr_result.c
RBR_API rbr_result * rbr_result_create(void);
RBR_API void rbr_result_destroy(rbr_result *out);

// rbr_param
// defined in rbr_param.c
RBR_API rbr_param * rbr_param_create(void);
RBR_API void rbr_param_destroy(rbr_param *par);

#ifdef __cplusplus
}
#endif

#endif
