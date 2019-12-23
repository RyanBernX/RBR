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
#include <getopt.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include "DCRBR.h"

#define BUFF_SIZE 1000
int cal_comm_num(int n, int k, int *labels);
void write_labels(const char *filename, int n, const int *labels);
void show_help();

int main(int argc, char **argv){
    adj_matrix A;
    DCRBR_out out;
    DCRBR_param param;
    int c, offset = 5;
    int k, comm_n;
    double tt = 0, obj_best = 0;
    char label_filename[BUFF_SIZE] = {0};
    int prob = 0, weighted = 0;

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
    param.maxIter = 25; param.verbose = 0; param.shuffle = 0; param.full = 0;
    param.funct_chk_freq = 10; param.tol = 4e-3;
    param.extract = RBR_ROUNDING;
    while ((c = getopt_long(argc, argv, "t:vhe:o:", long_opts, NULL)) != -1){
        switch (c){
            case 0:
                param.maxIter = atoi(optarg);
                break;
            case 1:
                offset = atoi(optarg);
                break;
            case 2:
                param.shuffle = atoi(optarg);
                break;
            case 3:
                param.full = 1;
                break;
            case 4:
                if (strcmp(optarg, "community") == 0){
                    prob = 0;
                } else if (strcmp(optarg, "maxcut") == 0){
                    prob = 1;
                } else {
                    fprintf(stderr, "Unrecognized problem: %s\n", optarg);
                    fprintf(stderr, "Valid value: community, maxcut.\n");
                    return 1;
                }
                break;
            case 5:
                weighted = 1;
                break;
            case 6:
                param.funct_chk_freq = atoi(optarg);
                break;
            case 't':
                param.tol = atof(optarg);
                break;
            case 'v':
                param.verbose = 1;
                param.k_param.verbose = 1;
                break;
            case 'h':
                show_help();
                return 0;
            case 'e':
                if (strcmp(optarg, "rounding") == 0){
                    param.extract = RBR_ROUNDING;
                } else if (strcmp(optarg, "kmeans") == 0){
                    param.extract = RBR_KMEANS;
                } else {
                    fprintf(stderr, "Unrecognized extraction method: %s\n", optarg);
                    fprintf(stderr, "Valid value: rounding, kmeans.\n");
                    return 1;
                } 
                break;
            case 'o':
                sprintf(label_filename, "%s", optarg);
                break;
        }
    }
    param.roundIter = param.maxIter - offset;

    if (argc - optind < 3 - param.full){
        show_help();
        return 1;
    }

    if (read_adj_matrix_csr(argv[optind], &A, 0, weighted))
        return 1;


    k = strtol(argv[optind+1], NULL, 10);
    if (param.full)
        param.p = k;
    else
        param.p = strtol(argv[optind+2], NULL, 10);

    srand((unsigned)time(NULL));

    switch(prob){
        case 0: /* community */
            for (int i = 0; i < 1; ++i){
                out = rbr(&A, k, param);
                if (i == 0 || out.Q > obj_best){
                    comm_n = cal_comm_num(A.n, k, out.labels);
                    if (label_filename[0] != '\0'){
                        write_labels(label_filename, A.n, out.labels);
                    }
                    obj_best = out.Q;
                }
                tt += out.elapsed;
                DCRBR_out_destroy(&out);
            }
            printf("nC: %d    Time: %f    CC: %f    S: %f    Q: %f\n", comm_n, tt, out.CC, out.S, obj_best);
            break;
        case 1: /* maxcut */
            out = rbr_maxcut(&A, k, param);
            printf("n: %d, iter: %d, cut: %e, time: %f\n", A.n, out.iter, out.funct_V, out.elapsed);
            DCRBR_out_destroy(&out);
            break;
        default:
            break;
    }

    adj_matrix_destroy(&A);
    return 0;
}

int cal_comm_num(int n, int k, int *labels){
    int *cnt = (int *)calloc(k, sizeof(int));
    int ret = 0;

    for (int i = 0; i < n; ++i)
        ++cnt[labels[i]];

    for (int i = 0; i < k; ++i)
        if (cnt[i] > 0) ++ret;

    free(cnt);
    return ret;
}

void write_labels(const char *filename, int n, const int *labels){
    FILE *fp;

    fp = fopen(filename, "w");

    if (fp == NULL){
        fprintf(stderr, "Cannot open file %s\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i)
        fprintf(fp, "%d %d\n", i, labels[i]);

    fclose(fp);
}


void show_help(){
    fprintf(stderr, "RBR help page\n");
    fprintf(stderr, "-----------------------------------------------------\n");
    fprintf(stderr, "[SYNOPSIS]\n");
    fprintf(stderr, "./rbr [OPTIONS] matrix k p\n");
    fprintf(stderr, "./rbr --full [OPTIONS] matrix k\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[ARGUMENTS]\n");
    fprintf(stderr, "matrix: Input matrix filename.\n");
    fprintf(stderr, "k: Number of clusters.\n");
    fprintf(stderr, "p: Number of non-zeros in U. Incompatible with --full\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]\n");
    fprintf(stderr, "-h,--help: Show help page.\n");
    fprintf(stderr, "-v,--verbose: Verbose output.\n");
    fprintf(stderr, "-e,--extract: Extraction method. Valid: rounding, kmeans. Default: rounding\n");
    fprintf(stderr, "--full: Use full storation of U.\n");
    fprintf(stderr, "--output file: Write the results into file.\n");
    fprintf(stderr, "--max-itr n: Max number of iterations. Default: 25.\n");
    fprintf(stderr, "--round-offset n: Perform rounding in last n iterations. Ignored when extraction method == kmeans. Default: 5.\n");
    fprintf(stderr, "--shuffle b: Optimize each row in a random order. Default: 0 (false).\n");
}
