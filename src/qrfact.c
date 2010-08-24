#include "qrfact.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, const char * argv[]){
	
  printf("----------Seccion 1 Numerov ----------\n");
  //double start = strtod(argv[1], (char **)NULL);
  double start=0.001;
  //double final = strtod(argv[2], (char **)NULL);
  double final=3.1416;
  int size = strtol(argv[1], (char **)NULL, 10);//atoi(argv[3]);//
	
  printf("Start:\t%lf\n",start);
  printf("Final:\t%lf\n", final);
  printf("Size:\t%d\n", size);
  printf("CPUs:\t%d\n",omp_get_num_procs());
  omp_set_num_threads(omp_get_num_procs()*2);

  double **T = NULL;
  printf("Matriz M\n");
  double **M = StartM(&T,size);
  PrintMatrix(M,size);

  printf("Matriz Q\n");
  double h;
  double **Q = StartQ(start,final,&h,size);
  PrintMatrix(Q,size);
  double **K = K_Process(M,Q,T,&h,size);
  printf("Matriz K\n");
  PrintMatrix(K,size);
  /* 
     M= I - (1/12)T
     K= (1/h^2)T + MQ
     Iterativo A = (M^(-1))*K
  */
  double **A= AllocateMatrixSpace(size);
  double **Qp= AllocateMatrixSpace(size);
  double **Mi= AllocateMatrixSpace(size);
  
  printf("\n");
  printf("Proceso Iterativo de QR\n");
  //PrintMatrix(A,size);
  //Qp = QR_Process(A);
	for (i=0;i<80;i++){
		A= Multiplication(M,K,size);
		Qp = QR_Method(A,size);
		K = K_Process(M,Qp,T,&h,size);
	}
		
  printf("----------Seccion 2 Corrimiento de Lambda------------\n");
	
  printf("----------Seccion 3 RK ----------\n");
	
	
}

double** AllocateMatrixSpace(int size){
  double **m = (double**)malloc(sizeof(double)*size);
  for(int i=0;i<size;i++)
    m[i]=(double*)malloc(sizeof(double)*size);
  return m;
}

double** StartM(double*** T, int size){
  double** Res = AllocateMatrixSpace(size);
  *T = AllocateMatrixSpace(size);
	
  for(int i=0;i<size;i++){
    (*T)[i][i]=2.0*(1.0/12.0);
    if(i>0)
      (*T)[i][i-1]=-1.0*(1.0/12.0);
    if(i<size-1)
      (*T)[i][i+1]=-1.0*(1.0/12.0);
  }

  //PrintMatrix(*T,size);

  for(int i=0;i<size;i++)
    for(int j=0;j<size;j++)
      Res[i][j] = (i==j ? 1.0 : 0.0) - (*T)[i][j];

  //PrintMatrix(Res,size);
  return Res; 
}  

double** StartQ(double Start, double End, double* h, int size){
  double** res = AllocateMatrixSpace(size);
  *h = (End-Start)/(double)size;
  double x= Start;

  for(int i=0; i<size; i++) {
    for (int j=0; j<size; j++){
      /* 
       * Evaluamos la función para generar Q
       * Monopolo Mágnetico: -1/(x^3)                                
       */
      if(i==j){
	res[i][j]=((-1.0)/(x*x*x));
	x=x+*h;
      }	else
	res[i][j]=0;
    }
  }
  
  return res;
  
}

double **K_Process(double **M, double **Q, double **T, double *h, int size) {
  // K= (1/h^2)*T+(MQ)
  double **MR = AllocateMatrixSpace(size);
  //double **P_MQ = AllocateMatrixSpace(size);
  //P_MQ = Multiplication(M,Q,size);
  double **P_MQ = Multiplication(M,Q,size);

  for(int i=0;i<size;i++)
    for(int j=0;j<size;j++)
      MR[i][j]= (((1/((*h)*(*h)))*T[i][j])+(P_MQ[i][j]));

  return MR;
}

double **Multiplication(double** matriz1, double** matriz2, int size){
  double **matrizR = AllocateMatrixSpace(size);

  for(int i=0;i<size;i++)
    for(int j=0;j<size;j++)
      for(int k=0;k<size;k++)
	matrizR[i][j]+=matriz1[i][k]*matriz2[k][j];

  //PrintMatrix(matrizR,size);
  return matrizR;
}

void PrintMatrix(double** matriz, int size) {
  for(int i=0;i<size;i++) {
    for(int j=0;j<size;j++)
      printf("%lf\t",matriz[i][j]);
    printf("\n");
  }
}



void Transpose(double **matriz){
	double tem;
	int i,j;
	for(i=0;i<SIZE;i++){
		for(j=i;j<SIZE;j++){
			if (i!=j) {
				tem=matriz[i][j];
				matriz[i][j]=matriz[j][i];
				matriz[j][i]=tem;
			}
		}
	}
}

double** QR_Method(double **a,int size){
	int i,j;
	double **q= AllocateMatrixSpace(size);
	double *columna_tem,tem;
	Transpose(a);
	CopyMatrix(a,q);
#pragma omp parallel for private(j)
	for(i=0;i<SIZE;i++){
		for(j=0;j<i;j++){
			tem=DotProduct(a[i],q[j]);
			tem/=DotProduct(q[j],q[j]);	
			columna_tem=DotProductVector(-tem,q[j]);
			columna_tem=SumVectors(q[i],columna_tem);
			CopyColumn(q[i],columna_tem);
		}
	}
	
#pragma omp parallel for private(tem)
	for(i=0;i<SIZE;i++) {
		tem=DotProduct(q[i],q[i]);
		DotProductVector_E(1.0/sqrt(tem),q[i]);
	}
	
#pragma omp parallel for private(j)
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			if (j>i) {
				r[i][j]=0;
			}else
				r[i][j]=DotProduct(a[i],q[j]);
		}
	}
	Transpose(a);
	Transpose(q);
	Transpose(r);
}

void CopyMatrix(double **a,double **b){
	int i,j;
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			b[i][j]=a[i][j];
		}
	}
}

double DotProduct(double *a,double *b){
	int i;
	double producto=0;
#pragma omp parallel for reduction(+:producto)
	for(i=0;i<SIZE;i++){
		producto+=a[i]*b[i];
	}
	return producto;
}

double *DotProductVector(double e,double *v){
	double *v_tem;
	int i;
	v_tem=(double*)malloc(sizeof(double)*SIZE);
#pragma omp parallel for 
	for(i=0;i<SIZE;i++){
		v_tem[i]=v[i]*e;
	}
	return v_tem;
}

void DotProductVector_E(double e,double *v){
	int i;
#pragma omp parallel for 
	for(i=0;i<SIZE;i++){
		v[i]=v[i]*e;
	}
}

double *SumVectors(double *a,double *b){
	int i;
	double *suma=(double*)malloc(sizeof(double)*SIZE);
#pragma omp parallel for 
	for(i=0;i<SIZE;i++){
		suma[i]=a[i]+b[i];
	}
	return suma;
}


void CopyColumn(double *a,double *b){
	int i;
#pragma omp parallel for 
	for(i=0;i<SIZE;i++){
		a[i]=b[i];
	}
}



void PrintMatrix(double **matriz){
	int i,j;
	printf("++++++++++++++++++++++++++++\n");
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			printf("%lf\t",matriz[i][j]);
		}
		printf("\n");
	}
	printf("++++++++++++++++++++++++++++\n");
}


void WriteToFile(FILE *fs, double** matriz, int n_col){
	int i,j;
	printf("Escritura de matriz en archivo\n");
	for(i=0;i<n_col;i++) {
		for(j=0;j<n_col;j++){
			fprintf(fs,"%lf\t",matriz[i][j]);
			printf("%lf ", matriz[i][j]);
		}
		fprintf(fs,"\n");
		printf("\n");
	}
	fclose(fs); 
}

*/
