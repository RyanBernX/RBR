/*
 * ===========================================================================
 *
 *       Filename:  rbr_result.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/28/2023 11:18:34 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#include "types.h"
#include "RBR.h"
#include <stdlib.h>

rbr_result * rbr_result_create(void){
    rbr_result * res = calloc(1, sizeof(rbr_result));
    return res;
}

void rbr_result_destroy(rbr_result *out){
    if (out == NULL) return;
    if (out->labels) free(out->labels);
    if (out->U)  free(out->U);
    if (out->iU) free(out->iU);

    free(out);
}

