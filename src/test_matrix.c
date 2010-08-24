#include "matrix.h"
#include <malloc.h>
#include <stdio.h>

int main(int argc, char** argv) {
  // Print some square matrix
  double** mat = (double**)calloc(3,sizeof(double*));
  for(int i=0; i<3; i++) {
    mat[i] = (double*)calloc(3,sizeof(double));
    for(int j=0; j<3; j++) {
      mat[i][j] = j+1;
    }
  }

  printf("Matrix\n");
  printSqrMat(mat,3);

  freeSqrMat(mat,3);


  // Create a Zero matrix
  mat = getZeroSqrMat(3);
  printf("Matrix\n");
  printSqrMat(mat,3);
  freeSqrMat(mat,3);


  // Compare two matrices
  mat = getZeroSqrMat(3);
  printf("Matrix1\n");
  printSqrMat(mat,3);

  double** mat1 = getZeroSqrMat(3);
  printf("Matrix2\n");
  printSqrMat(mat1,3);

  int eq = areEqualsSqrMat(mat,mat1,3);
  printf("Equals\t%d\n", eq);


  // Create an identity matrix
  mat = getIdentitySqrMat(3);
  printf("Matrix\n");
  printSqrMat(mat,3);
  freeSqrMat(mat,3);
  

  return 0;

}
