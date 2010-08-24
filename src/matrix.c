#include <malloc.h>

void applyTranspose(double** m, int n) {
   double tmp;

   for(int i=1; i<n; i++) {
      for(int j=0; j<i; j++) {
         tmp = m[i][j];
         m[i][j] = m[j][i];
         m[j][i] = tmp;
      }
   }
}

void printSqrMat(double** m, int size) {
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++)
      printf("%lf\t",m[i][j]);
    printf("\n");
  }
}
