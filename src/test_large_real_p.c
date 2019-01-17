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

int main(int argc, char **argv){
    adj_matrix A;
    DCRBR_out out;
    DCRBR_param param;
    int c, offset = 5;
    int k, comm_n;
    int threads;
    double tt = 0, obj_best = 0;
    char label_filename[BUFF_SIZE] = {0};

    /* parsing command line */
    struct option long_opts[] = {
        {"max-itr", 1, 0, 0},
        {"round-offset", 1, 0, 1},
        {"shuffle", 1, 0, 2},
        {"verbose", 0, 0, 'v'},
        {"extract", 1, 0, 'e'},
        {"output", 1, 0, 'o'}
    };
    param.maxIter = 25; param.verbose = 0; param.shuffle = 0;
    while ((c = getopt_long(argc, argv, "ve:o:", long_opts, NULL)) != -1){
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
            case 'v':
                param.verbose = 1;
                param.k_param.verbose = 1;
                break;
            case 'e':
                if (strcmp(optarg, "rounding") == 0){
                    param.extract = RBR_ROUNDING;
                } else if (strcmp(optarg, "kmeans") == 0){
                    param.extract = RBR_KMEANS;
                } else {
                    fprintf(stderr, "Unrecognized extraction method: %s\n", optarg);
                    fprintf(stderr, "Legal value: rounding, kmeans.\n");
                    return 1;
                } 
                break;
            case 'o':
                sprintf(label_filename, "%s", optarg);
                break;
        }
    }
    param.roundIter = param.maxIter - offset;

    if (argc - optind < 4){
        fprintf(stderr, "Usage ./main [OPTIONS] <csr_matrix> <k> <p> <threads>\n");
        return 1;
    }

    if (read_adj_matrix_csr(argv[optind], &A, 0))
        return 1;


    k = strtol(argv[optind+1], NULL, 10);
    param.p = strtol(argv[optind+2], NULL, 10);
    threads = strtol(argv[optind+3], NULL, 10);

    omp_set_num_threads(threads);

    srand((unsigned)time(NULL));
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

