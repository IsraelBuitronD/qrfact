#include "magnetic_monopole.h"
#include "vector.h"
#include "square_matrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

double** mm_numerov(double start, 
		    double final, 
		    int qr_iterations,
		    int matsize) {
  double **T = fillT(matsize);
  double **M = startM(T,matsize);
  double h = (final-start)/matsize;
  double **Q = startQ(start,final,matsize);
  double **K = kProcess(M,Q,T,h,matsize);
  double **Mi = getSqrMatInverse(M,matsize);
  double **Qp = NULL;
  double **A  = getSqrMatMulti(Mi,K,matsize);
  freeSqrMat(K,matsize);

  for(int i=0; i<qr_iterations; i++) {
    Qp = qrMethod(A,matsize);
    K  = kProcess(Mi,Qp,T,h,matsize);
    freeSqrMat(A,matsize);
    A  = getSqrMatMulti(Mi,K,matsize);

    freeSqrMat(Qp,matsize);
    freeSqrMat(K,matsize);
  }

  //freeSqrMat(Qp,matsize);
  freeSqrMat(Mi,matsize);
  //freeSqrMat(K,matsize);
  freeSqrMat(Q,matsize);
  freeSqrMat(M,matsize);
  freeSqrMat(T,matsize);

  return A;
}

double mm_lambda(double start,
		 double final,
		 int    size,
		 double delta,
		 double miu) {
  double miuf  = miu;
  int sigma_aux= 0;
  int sigma    = 0;
  double **P = NULL;
  double **T = fillT(size);
  double **M = startM(T,size);
  double   h = (final-start) / size;
  double **Q = startQ(start,final,size);
  double **K = kProcess(M,Q,T,h,size);
  //P= k-miu*M

  while(miu < 3.1416) {
    //printf("miu\t%lf\t",miu);
    P = pProcess(K,miu,M,size);
    sigma_aux = luMethod(P,size);
    freeSqrMat(P,size);
    //printf("sigaux\t%d\t",sigma_aux);
    //printf("sigma\t%d\n",sigma);
    
    if(sigma_aux > sigma){
      sigma = sigma_aux;
      miuf = miu;
      printf("#Valores propios\t%d\t", sigma);
      printf("Miu\t%lf\n",miuf);
    }
    miu += delta;
  }
  freeSqrMat(M,size);
  freeSqrMat(Q,size);
  freeSqrMat(T,size);
  freeSqrMat(K,size);

  return 0;
}

double mm_rungekutta(double start,
		     double final,
		     int    size,
		     double delta) {
  double k1,k2,k3,k4,l1,l2,l3,l4,lambda;
  double h = (final-start)/size;

  for(lambda = 0.001; lambda < 3.1416; lambda+=delta) {
    double u=0.0;
    double v=1.0;

    for(double x=0.001; x<final; x+=h) {
      k1=h*v;
      l1=h*function_uvx(u,x,lambda);
      k2=h*(v+(l1/2.0));
      l2=h*function_uvx((u+(k1/2.0)),(x+(h/2.0)),lambda);
      k3=h*(v+(l2/2.0));
      l3=h*function_uvx((u+(k2/2.0)),(x+(h/2.0)),lambda);
      k4=h*(v+(l3/2.0));
      l4=h*function_uvx((u+k3),(x+h),lambda);
      u+=((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));
      v+=((1.0/6.0)*(l1+(2.0*l2)+(2.0*l3)+l4));

      if(x==10.0)
	exit(0);
    }

    if(fabs(u)<0.05)
      printf("x:\t%lf\tlambda:\t%lf\n",fabs(u),lambda);

  }

  return lambda;
}

double function_uvx(double u, double x, double lambda) {
  return (((-1.0)/(x*x*x))-lambda)*u;
}

double** startQ(double start, double end, int size) {
  double** res = getZeroSqrMat(size);
  double h = (end-start)/size;
  double x = start;
 
  for(int i=0; i<size; i++) {
    for(int j=0; j<size; j++){
      /* 
       * Evaluate function to compute Q matrix.
       * Magnetic monopole function:
       * x^(-3)
       */
      if(i==j){
	res[i][j]=((-1.0)/(x*x*x));
	x+=h;
      }	else
	res[i][j]=0;
    }
  }
  
  return res;
}

double** kProcess(double **m, 
		  double **Q, 
		  double **T, 
		  double h, 
		  int size) {
  // K = (1/h^2)*T+(MQ)
  double **mr_matrix = getZeroSqrMat(size);
  double **mq_matrix = getSqrMatMulti(m,Q,size);
  int i,j;
  
#pragma omp parallel for private(j)
  for(i=0;i<size;i++)
    for( j=0;j<size;j++)
      mr_matrix[i][j] = ((1/(h*h)) * T[i][j] + 
			 mq_matrix[i][j]);

  freeSqrMat(mq_matrix, size);

  return mr_matrix;
}

double** qrMethod(double **a, int size) {
  int i,j;
  double *columna_tem,tem;
  applyTranspose(a,size);
  double **q = cloneSqrMat(a,size);

#pragma omp parallel for private(j)
  for(i=0;i<size;i++){
    for(j=0;j<i;j++){
      tem=DotProduct(a[i],q[j],size);
      tem/=DotProduct(q[j],q[j],size);	
      columna_tem=DotProductVector(-tem,q[j],size);
      columna_tem=SumVectors(q[i],columna_tem,size);
      CopyColumn(q[i],columna_tem,size);
    }
  }
	
#pragma omp parallel for private(tem)
  for(i=0;i<size;i++) {
    tem=DotProduct(q[i],q[i],size);
    applyScalarProduct(q[i],size,1.0/sqrt(tem));
  }
	
  double **r = getZeroSqrMat(size);
#pragma omp parallel for private(j)
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      r[i][j] = j>i ? 0 : DotProduct(a[i],q[j],size);

  applyTranspose(a,size);
  applyTranspose(q,size);
  applyTranspose(r,size);

  return q;
}


int luMethod(double **A, int size) {
  int k, sigma=0;
  double* U = (double*)calloc(size,sizeof(double));
  double* L = (double*)calloc(size,sizeof(double));
  double* a = (double*)calloc(size,sizeof(double));
  double* b = (double*)calloc(size,sizeof(double));
  double* c = (double*)calloc(size,sizeof(double));

  for(k=1; k<size; k++)
    a[k] = A[k-1][k];

  for(k=0; k<size; k++)
    b[k] = A[k][k];
  
  for(k=0; k<size-1; k++)
    c[k] = A[k+1][k];
  
  U[0] = A[0][0];
  
  if(U[0]<0)
    sigma++;
  
  for(k=0; k<size-1; k++) {
    L[k+1] = a[k+1] / U[k];
    U[k+1] = b[k+1] - L[k+1] * c[k];

    if(U[k+1]<0)
      sigma++;
  }
  
  free(U);
  free(L);
  free(a);
  free(b);
  free(c);
  
  return sigma;
}

double** pProcess(double** K, double miu, double** M, int size) {

  double **res = getZeroSqrMat(size);
  int i,j;

#pragma omp parallel for private(j)
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      res[i][j] = K[i][j] - miu * M[i][j];

  return res;
}

double** fillT(int size) {
  double** T = getZeroSqrMat(size);
  int i;

#pragma omp parallel for
  for(i=0; i<size; i++) {
    T[i][i]=2.0*(1.0/12.0);
    if(i>0)
      T[i][i-1]=-1.0*(1.0/12.0);
    if(i<size-1)
      T[i][i+1]=-1.0*(1.0/12.0);
  }

  return T;
}

double** startM(double** T, int size) {
  double** res = getZeroSqrMat(size);

  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      res[i][j] = (i==j ? 1.0 : 0.0) - T[i][j];

  return res;
}
