void kb07ad_(double*, int*, int*);


void dsort(int n, double *data, int *index){
  kb07ad_(data, &n, index);
  for (int i = 0; i < n; ++i) --index[i];
}
