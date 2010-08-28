#include "square_matrix.h"
#include <stdlib.h>
#include <math.h>
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
  

  // Compute matrix determinant 
  mat = getZeroSqrMat(4);
  mat[0][0]=1;
  mat[0][1]=2;
  mat[0][2]=3;
  mat[0][3]=4;
  mat[1][0]=4;
  mat[1][1]=1;
  mat[1][2]=2;
  mat[1][3]=3;
  mat[2][0]=3;
  mat[2][1]=4;
  mat[2][2]=1;
  mat[2][3]=2;
  mat[3][0]=2;
  mat[3][1]=3;
  mat[3][2]=4;
  mat[3][3]=1;
  printf("Matrix\n");
  printSqrMat(mat,4);
  double det = determinant(mat, 4);
  printf("Det\t%lf\t%s\n", det, (det==-160.0 ? "correct" : "wrong") );
  freeSqrMat(mat,4);


  // Compute cofactor matrix
  /* mat = getZeroSqrMat(3); */
  /* mat[0][0]=-4; */
  /* mat[0][1]=2; */
  /* mat[0][2]=1; */
  /* mat[1][0]=1; */
  /* mat[1][1]=-10; */
  /* mat[1][2]=1; */
  /* mat[2][0]=3; */
  /* mat[2][1]=-1; */
  /* mat[2][2]=-2; */
  /* double** comat = comatrix(mat,3); */
  /* printf("CoMatrix\n"); */
  /* printSqrMat(comat,3); */
  /* freeSqrMat(comat,3); */

  mat = getZeroSqrMat(3);
  mat[0][0]=1;
  mat[0][1]=2;
  mat[0][2]=0;
  mat[1][0]=-1;
  mat[1][1]=1;
  mat[1][2]=1;
  mat[2][0]=1;
  mat[2][1]=2;
  mat[2][2]=3;

  printf("Matrix\n");
  printSqrMat(mat,3);
  
  det = determinant(mat, 3);
  printf("Deter\t%lf\n", det );

  double** comat = comatrix(mat,3);
  printf("CoMatrix\n");
  printSqrMat(comat,3);

  applyTranspose(comat,3);
  printf("Apply Transpose CoMatrix\n");
  printSqrMat(comat,3);

  applyScalarMultiplication(comat,3,1/det);
  printf("Apply Scalar Multi\n");
  printSqrMat(comat,3);

  freeSqrMat(comat,3);




  return 0;

}

int testSqrMatInverse() {
  // Compute inverse matrix
  int size = 3;
  double **mat = getZeroSqrMat(size);
  mat[0][0]=2;
  mat[0][1]=0;
  mat[0][2]=1;
  mat[1][0]=3;
  mat[1][1]=0;
  mat[1][2]=0;
  mat[2][0]=5;
  mat[2][1]=1;
  mat[2][2]=1;
  
  printf("Matrix\n");
  printSqrMat(mat,size);

  double **comat = getSqrMatInverse(mat,size);
  printf("Inverse Matrix\n");
  printSqrMat(comat,size);

  freeSqrMat(mat,size);
  freeSqrMat(comat,size);

  return 0;
}
