#include <stdio.h>

void dsort(int, double*, int*);

int main(int argc, char **argv){
  double a[] = {0.4, 0.1, 0.3, 0.5, 0.2};
  int ia[5];

  dsort(5, a, ia);
  for (int i = 0; i < 5; ++i)
    printf("(%.2f, %d)\n", a[i], ia[i]);

  return 0;
}