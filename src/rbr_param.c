/*
 * ===========================================================================
 *
 *       Filename:  rbr_param.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/29/2023 11:46:38 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#include "RBR.h"
#include "types.h"
#include <stdlib.h>

rbr_param * rbr_param_create(void){
    rbr_param * param = (rbr_param*)calloc(1, sizeof(rbr_param));

    if (param == NULL){
        return NULL;
    }

    param->maxIter = 25;
    param->verbose = 0;
    param->shuffle = 0;
    param->full = 0;
    param->funct_chk_freq = 10;
    param->tol = 4e-3;
    param->extract = RBR_ROUNDING;

    param->k_param.verbose = 0;
    param->k_param.init = KMEANS_SAMPLE;
    param->k_param.maxit = 20;

    return param;
}

void rbr_param_destroy(rbr_param *p){
    if (p != NULL){
        free(p);
    }
}

