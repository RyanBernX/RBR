#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <DCRBR.h>


void dsort(int n, double *data, int *index);
double quickselect(double *arr, int n, int k, int *index);

int main(int argc, char **argv){
  int n, k;

  n = strtol(argv[1], NULL, 10);
  k = strtol(argv[2], NULL, 10);

  int index[n], iPos = 0;
  //double data[n], data_t[n];

  double data[] = {-5, 4, -4, -1, 2, 0};
  solve_sub_U(n, k, data, index, &iPos);

  printf("%d\n", iPos);
  for (int i = 0; i < iPos; ++i)
    printf("(%d %f)\n", index[i], data[i]);

/*
  srand((unsigned)time(0));
  for (int i = 0; i < n; ++i)
    data[i] = 1.0 * rand() / RAND_MAX;

  t1 = clock();
  for (int i = 0; i < 10000; ++i){
    memcpy(data_t, data, n * sizeof(double));
    quickselect(data_t, n, k, index);
  }
  t1 = (clock() - t1) / CLOCKS_PER_SEC;

  t2 = clock();
  for (int i = 0; i < 10000; ++i){
    memcpy(data_t, data, n * sizeof(double));
    dsort(n, data, index);
  }
  t2 = (clock() - t2) / CLOCKS_PER_SEC;

  t3 = clock();
  for (int i = 0; i < 10000; ++i){
    memcpy(data_t, data, n * sizeof(double));
    nth_element(data_t, k, n, index);
  }
  t3 = (clock() - t3) / CLOCKS_PER_SEC;

  printf("METHOD    TIME\n");
  printf("sort      %f\n", t2);
  printf("select    %f\n", t1);
  printf("nth       %f\n", t3);
*/
  return 0;

}
