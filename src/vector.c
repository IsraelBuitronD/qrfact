#include "vector.h"
#include "square_matrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

double DotProduct(double* a, double* b, int size) {
  int i;
  double producto=0;
  
#pragma omp parallel for reduction(+:producto)
  for(i=0;i<size;i++)
    producto+=a[i]*b[i];

  return producto;
}

double *DotProductVector(double e, double *v, int size) {
  int i;
  double *v_tem = (double*)calloc(size, sizeof(double));

#pragma omp parallel for 
  for(i=0; i<size; i++)
    v_tem[i]=v[i]*e;

  return v_tem;
}

double *SumVectors(double *a,double *b,int size) {
  int i;
  double *suma=(double*)calloc(size, sizeof(double));

#pragma omp parallel for 
  for(i=0; i<size; i++)
    suma[i]=a[i]+b[i];
  
  return suma;
}

void CopyColumn(double* src, double* des, int size) {
  int i;

#pragma omp parallel for 
  for(i=0; i<size; i++)
    src[i] = des[i];

}

double* applyScalarProduct(double *v, int size, double e) {
  int i;

#pragma omp parallel for 
  for(i=0; i<size; i++)
    v[i]=v[i]*e;

  return v;
}
 
