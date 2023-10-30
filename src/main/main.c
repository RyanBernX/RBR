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

#include "RBR.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <string.h>

#define BUFF_SIZE 1000
RBR_INT cal_comm_num(RBR_INT n, RBR_INT k, const RBR_INT *labels);
void write_labels(const char *filename, RBR_INT n, const RBR_INT *labels);
void show_help(void);

int main(int argc, char **argv){
    int c, offset = 5;
    char label_filename[BUFF_SIZE] = {0};
    int prob = 0, weighted = 0;

    rbr_param *param = rbr_param_create();

    if (param == NULL){
        fprintf(stderr, "[ERROR] rbr: failed to create parameter struct\n");
        return 1;
    }

    /* parsing command line */
    struct option long_opts[] = {
        {"max-itr", 1, 0, 0},
        {"round-offset", 1, 0, 1},
        {"shuffle", 1, 0, 2},
        {"full", 0, 0, 3},
        {"prob", 1, 0, 4},
        {"weighted", 0, 0, 5},
        {"funct-chk-freq", 1, 0, 6},
        {"tol", 1, 0, 't'},
        {"verbose", 0, 0, 'v'},
        {"help", 0, 0, 'h'},
        {"extract", 1, 0, 'e'},
        {"output", 1, 0, 'o'}
    };

    while ((c = getopt_long(argc, argv, "t:vhe:o:", long_opts, NULL)) != -1){
        switch (c){
            case 0:
                param->maxIter = atoi(optarg);
                break;
            case 1:
                offset = atoi(optarg);
                break;
            case 2:
                param->shuffle = atoi(optarg);
                break;
            case 3:
                param->full = 1;
                break;
            case 4:
                if (strcmp(optarg, "community") == 0){
                    prob = 0;
                } else if (strcmp(optarg, "maxcut") == 0){
                    prob = 1;
                } else {
                    fprintf(stderr, "Unrecognized problem: %s\n", optarg);
                    fprintf(stderr, "Valid value: community, maxcut.\n");
                    rbr_param_destroy(param);
                    return 1;
                }
                break;
            case 5:
                weighted = 1;
                break;
            case 6:
                param->funct_chk_freq = atoi(optarg);
                break;
            case 't':
                param->tol = atof(optarg);
                break;
            case 'v':
                param->verbose = 1;
                param->k_param.verbose = 1;
                break;
            case 'h':
                show_help();
                return 0;
            case 'e':
                if (strcmp(optarg, "rounding") == 0){
                    param->extract = RBR_ROUNDING;
                } else if (strcmp(optarg, "kmeans") == 0){
                    param->extract = RBR_KMEANS;
                } else {
                    fprintf(stderr, "Unrecognized extraction method: %s\n", optarg);
                    fprintf(stderr, "Valid value: rounding, kmeans.\n");
                    rbr_param_destroy(param);
                    return 1;
                } 
                break;
            case 'o':
                sprintf(label_filename, "%s", optarg);
                break;
        }
    }
    param->roundIter = param->maxIter - offset;

    if (argc - optind < 3 - param->full){
        fprintf(stderr, "Please specify the parameter p (sparsity argument).\n");
        return 1;
    }

    adj_matrix *A = read_adj_matrix_csr(argv[optind], 0, weighted);

    if (A == NULL){
        fprintf(stderr, "[ERROR] Unable to open file \'%s\'\n", argv[optind]);
        return 1;
    }


    RBR_INT k = strtol(argv[optind+1], NULL, 10);
    if (param->full){
        param->p = k;
    } else {
        param->p = strtol(argv[optind+2], NULL, 10);
    }

    srand((unsigned)time(NULL));

    rbr_result * out = NULL;
    switch(prob){
        case 0: /* community */
            out = rbr(A, k, param);
            RBR_INT comm_n = cal_comm_num(A->n, k, out->labels);
            if (label_filename[0] != '\0'){
                write_labels(label_filename, A->n, out->labels);
            }
            printf("nC: %d    Time: %f    CC: %f    S: %f    Q: %f\n", (int)comm_n, out->elapsed, out->CC, out->S, out->Q);
            break;
        case 1: /* maxcut */
            out = rbr_maxcut(A, k, param);
            printf("n: " _fmt_int "    iter: %d    cut: %e    time: %f\n", A->n, out->iter, out->funct_V, out->elapsed);
            break;
        default:
            break;
    }

    rbr_param_destroy(param);
    rbr_result_destroy(out);
    adj_matrix_destroy(A);
    return 0;
}

RBR_INT cal_comm_num(RBR_INT n, RBR_INT k, const RBR_INT *labels){
    RBR_INT *cnt = (RBR_INT*)calloc(k, sizeof(RBR_INT));
    RBR_INT ret = 0;

    if (cnt == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return -1;
    }

    for (RBR_INT i = 0; i < n; ++i)
        ++cnt[labels[i]];

    for (RBR_INT i = 0; i < k; ++i)
        if (cnt[i] > 0) ++ret;

    free(cnt);
    return ret;
}

void write_labels(const char *filename, RBR_INT n, const RBR_INT *labels){
    FILE *fp;

    fp = fopen(filename, "w");

    if (fp == NULL){
        fprintf(stderr, "Cannot open file %s\n", filename);
        return;
    }

    for (RBR_INT i = 0; i < n; ++i){
        fprintf(fp, _fmt_int " " _fmt_int "\n", i, labels[i]);
    }

    fclose(fp);
}


void show_help(void){
    fprintf(stderr, "RBR help page\n");
    fprintf(stderr, "-----------------------------------------------------\n");
    fprintf(stderr, "[SYNOPSIS]\n");
    fprintf(stderr, "    ./rbr [OPTIONS] matrix k p\n");
    fprintf(stderr, "    ./rbr --full [OPTIONS] matrix k\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[ARGUMENTS]\n");
    fprintf(stderr, "    matrix: Input matrix filename.\n");
    fprintf(stderr, "         k: Number of clusters.\n");
    fprintf(stderr, "         p: Number of non-zeros in U. Incompatible with --full\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]\n");
    fprintf(stderr, "    -h,--help: Show help page.\n");
    fprintf(stderr, "    -v,--verbose: Verbose output.\n");
    fprintf(stderr, "    -e,--extract: Extraction method. Valid: rounding, kmeans. Default: rounding\n");
    fprintf(stderr, "    --full: Use full storation of U.\n");
    fprintf(stderr, "    --output file: Write the results into file.\n");
    fprintf(stderr, "    --max-itr n: Max number of iterations. Default: 25.\n");
    fprintf(stderr, "    --round-offset n: Perform rounding in last n iterations. Ignored when extraction method == kmeans. Default: 5.\n");
    fprintf(stderr, "    --shuffle b: Optimize each row in a random order. Default: 0 (false).\n");
}

