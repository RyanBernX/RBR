#include <algorithm>
#include <numeric>
#include <string.h>
#include "DCRBR.h"

extern "C" {

void solve_sub_U(int n, int p, double *a, int *ia, int *iPos) {
  //for (int i = 0; i < n; ++i) index[i] = i;
  std::iota(ia, ia + n, 0);

  /* find the position of nth_element */  
  std::nth_element(ia, ia + p - 1, ia + n, [&a](int i1, int i2) { return *(a + i1) < *(a + i2); } );

  /* discard those positive or zero entries */
  double *a_tmp = (double*)calloc(p, sizeof(double));
  int *ia_tmp = (int*)calloc(p, sizeof(int));
  double min_pos;
  int imin_pos = -1, num = 0;
  for (int i = 0; i < p; ++i){
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
    memcpy(a, a_tmp, num * sizeof(double));
    memcpy(ia, ia_tmp, num * sizeof(int));
  }

  free(a_tmp); free(ia_tmp);

  *iPos = num;
}

}

