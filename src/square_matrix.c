#include "square_matrix.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

void applyTranspose(double** m, int n) {
  double tmp;
  int i,j;

#pragma omp parallel for private(j,tmp)
  for(i=1; i<n; i++) {
    for(j=0; j<i; j++) {
      tmp = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = tmp;
    }
  }
}

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
  int i,j;

  for(i=0; i<size; i++)
    for(j=0; j<size; j++)
      if(m1[i][j] != m2[i][j])
	return NOT_EQUAL_VALUES;

  return EQUAL_VALUES;
}

double** getZeroSqrMat(int size) {
  double** mat = (double**)calloc(size,sizeof(double*));

  int i;
#pragma omp parallel for 
  for(i=0; i<size; i++)
    mat[i] = (double*)calloc(size,sizeof(double));

  return mat;
}

double** freeSqrMat(double** m, int size) {
  int i;
#pragma omp parallel for 
  for(i=0; i<size; i++)
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

double** comatrix(double **a, int size) {
  double **b = getZeroSqrMat(size);
  double **c = getZeroSqrMat(size-1);

  for(int j=0;j<size;j++) {
    for(int i=0;i<size;i++) {

      // Form the adjoint a_ij
      for(int ii=0, i1 = 0; ii<size; ii++) {
	if (ii == i)
	  continue;
	for(int jj=0, j1=0; jj<size; jj++) {
	  if (jj == j)
	    continue;
	  c[i1][j1] = a[ii][jj];
	  j1++;
	}
	i1++;
      }

      // Compute determinat
      double det = determinant(c,size-1);

      // Fill comatrix elements
      b[i][j] = pow(-1.0,i+j+2.0) * det;
    }
  }

  freeSqrMat(c,size-1);

  return b;
}

double** applyScalarMultiplication(double** m, int size, double scalar) {
  int i,j;
  
#pragma omp parallel for private(j)
  for(i=0; i<size; i++)
    for(j=0; j<size; j++)
      m[i][j] *= scalar;
  
  return m;
}

double** getSqrMatInverse(double** m, int size) {
  double det = determinant(m, 3);

  if(det==0.0) {
    fprintf(stderr, "Determinant zero implies not inverse matrix.\n");
    exit(EXIT_FAILURE);
  }

  double** comat = comatrix(m,size);
  applyTranspose(comat,size);
  applyScalarMultiplication(comat,size,1/det);

  return comat;
}

double** cloneSqrMat(double **a, int size) {
  double **b = getZeroSqrMat(size);
  int i,j;
  
#pragma omp parallel for private(j)
  for(i=0; i<size; i++)
    for(j=0; j<size; j++)
      b[i][j]=a[i][j];
    
  return b;
}
