#include "matrix.h"
#include <stdlib.h>

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

/*
double** getTranspose(double** m, int n) {
  double** transpose = (double**)calloc(n,sizeof(double*));

  // TODO: implement this function
  
  return transpose;
}
*/

void printSqrMat(double** m, int size) {
  printDelimSqrMat(m,'\t',size);
}

void printDelimSqrMat(double** m, char delim, int size) {
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++)
      printf("%lf%c",m[i][j],delim);
    printf("\n");
  }
}

int areEqualsSqrMat(double** m1, double** m2, int n) {
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(m1[i][j] != m2[i][j])
	return NOT_EQUAL_VALUES;

  return EQUAL_VALUES;
}

double** getZeroSqrMat(int size) {
  double** mat = (double**)calloc(3,sizeof(double*));
  for(int i=0; i<3; i++)
    mat[i] = (double*)calloc(3,sizeof(double));

  return mat;
}

double** freeSqrMat(double** m, int size) {
  for(int i=0; i<3; i++)
    free(m[i]);
  free(m);
  return m;
}

double** getIdentitySqrMat(int size) {
  double** mat = (double**)calloc(3,sizeof(double*));
  for(int i=0; i<3; i++) {
    mat[i] = (double*)calloc(3,sizeof(double));
    mat[i][i] = 1;
  }

  return mat;
}
