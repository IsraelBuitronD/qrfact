#include "qrfact.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

double global[SIZE][SIZE]={{12,-51,4},{6,167,-68},{-4,24,-41}};

double **q,**r;

int main(){
  int i,hilos=0;
  double **a,**a_tem;

  printf("CPUs:\t%d\n",omp_get_num_procs());
  hilos=omp_get_num_procs()*2;
  omp_set_num_threads(4);
  printf("hilos: %d\n",omp_get_num_threads());
  a_tem=pedirMemoriaMatriz();
  q=pedirMemoriaMatriz();r=pedirMemoriaMatriz();
  a=crearA(0.0,2.0,SIZE);
  //copiarMatriz(a,q);
  for (i=0;i<20;i++){
    qr(a);
    transpuesta(q);
    a=multiplica(q,a);
    transpuesta(q);
    a=multiplica(a,q);
  }
  printf("Imprimiendo Q:\n");
  imprime(a);
  printf("Imprimiendo Q:\n");
  imprime(q);
  printf("Imprimiendo R:\n");
  imprime(r);
}

void qr(double **a){
  int i,j;
  double *columna_tem,tem;
  transpuesta(a);
  copiarMatriz(a,q);
#pragma omp parallel for private(j)
  for(i=0;i<SIZE;i++){
    for(j=0;j<i;j++){
      tem=productoPunto(a[i],q[j]);
      tem/=productoPunto(q[j],q[j]);	
      columna_tem=productoEscalarVector(-tem,q[j]);
      columna_tem=sumaVectores(q[i],columna_tem);
      copiarColumna(q[i],columna_tem);
    }
  }
	
#pragma omp parallel for private(tem)
  for(i=0;i<SIZE;i++) {
    tem=productoPunto(q[i],q[i]);
    productoEscalarVector1(1.0/sqrt(tem),q[i]);
  }
	
#pragma omp parallel for private(j)
  for(i=0;i<SIZE;i++){
    for(j=0;j<SIZE;j++){
      r[i][j]=productoPunto(a[i],q[j]);
    }
  }
  transpuesta(a);
  transpuesta(q);
  transpuesta(r);
}

double **pedirMemoriaMatriz(){
  int i;
  double **m;
  m=(double**)malloc(sizeof(double)*SIZE);
  for(i=0;i<SIZE;i++){
    m[i]=(double*)malloc(sizeof(double)*SIZE);
  }
  return m;
}

double *pedirMemoriaVector(){
  double *v;
  v=(double*)malloc(sizeof(double)*SIZE);
  return v;
}

void inicializaA(double **a){
  int i,j;
  for(i=0;i<SIZE;i++){
    for(j=0;j<SIZE;j++){
      a[i][j]=global[i][j];
    }
  }
}

void copiarMatriz(double **a,double **b){
  int i,j;
  for(i=0;i<SIZE;i++){
    for(j=0;j<SIZE;j++){
      b[i][j]=a[i][j];
    }
  }
}

void transpuesta(double **matriz){
  double tem;
  int i,j;
#pragma parallel for private(j,tem)
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

void imprime(double **matriz){
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

double productoPunto(double *a,double *b){
  int i;
  double producto=0;
#pragma omp parallel for reduction(+:producto)
  for(i=0;i<SIZE;i++){
    producto+=a[i]*b[i];
  }
  return producto;
}

void productoEscalarVector1(double e,double *v){
  int i;
#pragma omp parallel for 
  for(i=0;i<SIZE;i++){
    v[i]=v[i]*e;
  }
}

double *productoEscalarVector(double e,double *v){
  double *v_tem;
  int i;
  v_tem=pedirMemoriaVector();
#pragma omp parallel for 
  for(i=0;i<SIZE;i++){
    v_tem[i]=v[i]*e;
  }
  return v_tem;
}

double *sumaVectores(double *a,double *b){
  int i;
  double *suma=pedirMemoriaVector();
#pragma omp parallel for 
  for(i=0;i<SIZE;i++){
    suma[i]=a[i]+b[i];
  }
  return suma;
}

void copiarColumna(double *a,double *b){
  int i;
#pragma omp parallel for 
  for(i=0;i<SIZE;i++){
    a[i]=b[i];
  }
}

double **crearA(double intervaloInicio, double intervaloFinal, double numParticiones){
  double **resultado;
  resultado=pedirMemoriaMatriz();
  int i,j;
  double tamParticion=(intervaloFinal-intervaloInicio)/numParticiones;
  double T,Q,div,actual=intervaloInicio;
  div=1/tamParticion*tamParticion;
  for(i=0;i<numParticiones;i++){
    for(j=0;j<numParticiones;j++){
      if(i==j){
	T=2;
	actual=intervaloInicio+tamParticion*i;
	Q=(actual-0.5)*(actual-0.5);
      }
      else{
	if( (i<numParticiones-1&&(i==j+1||i==j-1))||(i>0&&(i==j-1||i==j+1 ) )){
	  T=-1;
	}
	else
	  T=0;
	Q=0;
      }
      resultado[i][j]=div*T+Q;
    }
  }
  return resultado;
}

double **multiplica(double **matriz1,double **matriz2){
  int i,j,k;
  double **matrizR=pedirMemoriaMatriz();
  for(i=0;i<SIZE;i++){
    for(j=0;j<SIZE;j++){
      for(k=0;k<SIZE;k++){
	matrizR[i][j]+=matriz1[i][k]*matriz2[k][j];
      }
    }
  }
  return matrizR;
}

double **readFromFile(FILE *fe){
  double **A;
  int i,j,n_col=0;
  if (fe == NULL) {
    fprintf(stderr, "Can't open input file\n");
    exit(EXIT_FAILURE);
  } else {
    printf("Lectura de la matriz en archivo\n");
    fscanf(fe,"%d",&n_col);
    A=(double **)malloc (n_col*sizeof(double*));
    for(i=0;i<n_col;i++)
      A[i]=(double *)malloc (n_col*sizeof(double));
    for(i=0;i<n_col;i++) {
      for(j=0;j<n_col;j++){
	fscanf(fe, "%lf", &A[i][j]);
	printf("%lf ", A[i][j]);
      }
      printf("\n");
    }
  }
  fclose(fe);
  return A;	
}

void writeToFile(FILE *fs, double** matriz, int n_col){
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

