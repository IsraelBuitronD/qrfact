#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

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

int areEqualsSqrMat(double** m1, double** m2, int size) {
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      if(m1[i][j] != m2[i][j])
	return NOT_EQUAL_VALUES;

  return EQUAL_VALUES;
}

double** getZeroSqrMat(int size) {
  double** mat = (double**)calloc(size,sizeof(double*));
  for(int i=0; i<size; i++)
    mat[i] = (double*)calloc(size,sizeof(double));

  return mat;
}

double** freeSqrMat(double** m, int size) {
  for(int i=0; i<size; i++)
    free(m[i]);
  free(m);
  return m;
}

double** getIdentitySqrMat(int size) {
  double** mat = (double**)calloc(size,sizeof(double*));
  for(int i=0; i<size; i++) {
    mat[i] = (double*)calloc(size,sizeof(double));
    mat[i][i] = 1;
  }

  return mat;
}

double determinant(double **a,int size) {
  if (size < 1) { // Invalid matrix size
    fprintf(stderr, "Matrix size must be greater than zero.\n");
    exit(EXIT_FAILURE);
  } 
  
  if(size == 1) // 1x1 Matrix
    return a[0][0];
  
  if(size == 2) // 2x2 Matrix
    return a[0][0] * a[1][1] - a[1][0] * a[0][1];

  // nxn Matrix => n>=3
  
  double det = 0;
  
  for(int j1=0; j1<size; j1++) {
    double **m = getZeroSqrMat(size-1);

    for(int i=1;i<size;i++)
      for(int j=0, j2=0; j<size; j++) {
	if (j == j1)
	  continue;
	m[i-1][j2] = a[i][j];
	j2++;
      }

    det += pow(-1.0,j1+2.0) * a[0][j1] * determinant(m,size-1);

    freeSqrMat(m,size-1);    
  }

  return det;
}
