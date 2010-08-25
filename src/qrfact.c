#include "qrfact.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "matrix.h"
int main(int argc, const char * argv[]){
  int i;	
  printf("----------Seccion 1 Numerov ----------\n");
  //double start = strtod(argv[1], (char **)NULL);
  double start=0.001;
  //double final = strtod(argv[2], (char **)NULL);
  double final=3.1416;
  int size = strtol(argv[1], (char **)NULL, 10);//atoi(argv[3]);//
  printf("\tProgramas para evaluar el Monopolo Mágnetico: -1/(x^3)\n");	
  printf("Start:\t%lf\n",start);
  printf("Final:\t%lf\n", final);
  printf("Size:\t%d\n", size);
  printf("CPUs:\t%d\n",omp_get_num_procs());
  omp_set_num_threads(omp_get_num_procs()*2);

  double **T = NULL;
  printf("Matriz M\n");
  double **M = StartM(&T,size);
  //PrintMatrix(M,size);

  printf("Matriz Q\n");
  double h;
  double **Q = StartQ(start,final,&h,size);
  //PrintMatrix(Q,size);
  double **K = K_Process(M,Q,T,&h,size);
  printf("Matriz K\n");
  //PrintMatrix(K,size);
  /* 
     M= I - (1/12)T
     K= (1/h^2)T + MQ
     Iterativo A = (M^(-1))*K
  */
  double **A= AllocateMatrixSpace(size);
  double **Qp=AllocateMatrixSpace(size);
  double **Mi= AllocateMatrixSpace(size);
  Mi= getSqrMatInverse(M,size);
  double eps= 0.400;
 //printf("\t---Matriz Transpuesta M\n");
 //PrintMatrix(Mi,size);
  printf("\n");
  printf("Proceso Iterativo de QR\n");
  //PrintMatrix(A,size);
  //Qp = QR_Process(A);
    A= Multiplication(Mi,K,size);
	for(i=0;i<100;i++){//while(*A[1,1]>eps){ //Cambiar por un while
		Qp = QR_Method(A,size);
		//printf("-----Qp-----\n");
		//PrintMatrix(Qp,size);
		K = K_Process(Mi,Qp,T,&h,size);
		//printf("-----K-----\n");
		//PrintMatrix(K,size);
		A= Multiplication(Mi,K,size);
		//printf("-----A-----\n");
		//PrintMatrix(A,size);
		printf("Iteración %d",i);	
	}
		printf("-----A----- After Iteration\n");
  		PrintMatrix(A,size);
		freeSqrMat(A,size );
	    freeSqrMat(Mi,size );
		freeSqrMat(Qp,size );
	  printf("\n\n\n----------Seccion 2 Corrimiento de Lambda------------\n");
		double **P=AllocateMatrixSpace(size);
		double miu=-15.0;
		double miuf;
		double delta=0.0001;
		int sigma_aux=0;
		int sigma=0; 
		K = K_Process(M,Q,T,&h,size);
		//P= k-miu*M
		printf("Proceso LU\n");
		miuf=miu;
		while(miu<3.1416){
			P= P_Process(K,miu,M,size);
			sigma_aux=LU_Method(P,size);
			if(sigma_aux>sigma){
				sigma= sigma_aux;
				miuf=miu;
				printf("\n\tNo. de Valores propios son: %d  para un \n ", sigma);
				printf("Miu de: %lf",miuf);
			}
			miu+=delta;
		}
		freeSqrMat(M,size);
		freeSqrMat(Q,size);
		freeSqrMat(T,size);
		freeSqrMat(K,size);
		freeSqrMat(P,size);
  printf("\n\n\n----------Seccion 3 RK ----------\n");
		double lambda=0.001;
		double x,u,v,k1,k2,k3,k4,l1,l2,l3,l4;
		delta = 0.001;
		printf("\tH %lf",h);
		while (lambda<3.1416){
			u=0.0;v=1.0;
			for(x=0.001;x<final;x+=h)
			{
				k1=h*v;
				//printf("k1 %lf ",k1);
				//printf("Entrada de f: %lf  %lf  %lf\n",u,x,lambda);
				l1=h*fuvx(u,x,lambda);
			    //printf("l1 %lf ",l1);
				k2=h*(v+(l1/2.0));
				//printf("k2 %lf ",k2);
				l2=h*fuvx((u+(k1/2.0)),(x+(h/2.0)),lambda);
				//printf("l2 %lf ",l2);
				k3=h*(v+(l2/2.0));
				//printf("k3 %lf ",k3);
				l3=h*fuvx((u+(k2/2.0)),(x+(h/2.0)),lambda);
				//printf("l3 %lf ",l3);
				k4=h*(v+(l3/2.0));
				//printf("k4 %lf ",k4);
				l4=h*fuvx((u+k3),(x+h),lambda);
				//printf("l4 %lf ",l4);
				u=u+((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));
				//printf("u %lf ",u);
				v=v+((1.0/6.0)*(l1+(2.0*l2)+(2.0*l3)+l4));
				//printf("v %lf \n\n",v);
				if(x==10.0)
					exit(0);
			}
			//printf("\nLambda:%lf %lf",fabs(u),lambda);
			if(fabs(u)<0.05)
				printf("\nLambda:%lf %lf",fabs(u),lambda);
			lambda+=delta;
		}	
		//printf("\tendwhile\n");
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
#pragma omp parallel for	
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
int i,j;
#pragma omp parallel for private(j)
  for( i=0;i<size;i++)
    for( j=0;j<size;j++)
      MR[i][j]= (((1/((*h)*(*h)))*T[i][j])+(P_MQ[i][j]));

  return MR;
}

double **Multiplication(double** matriz1, double** matriz2, int size){
  double **matrizR = AllocateMatrixSpace(size);
int i,j,k;
//#pragma omp parallel for private(j,k)
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      for(k=0;k<size;k++)
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


double** QR_Method(double **a,int size){
	int i,j;
	double **q= AllocateMatrixSpace(size);
	double **r= AllocateMatrixSpace(size);
	double *columna_tem,tem;
	applyTranspose(a,size);
	CopyMatrix(a,q,size);
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
		DotProductVector_E(1.0/sqrt(tem),q[i],size);
	}
	
#pragma omp parallel for private(j)
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			if (j>i) {
				r[i][j]=0;
			}else
				r[i][j]=DotProduct(a[i],q[j],size);
		}
	}
	applyTranspose(a,size);
	applyTranspose(q,size);
	applyTranspose(r,size);
	return q;
}

void CopyMatrix(double **a,double **b,int size){
	int i,j;
	#pragma omp parallel for private(j)
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			b[i][j]=a[i][j];
		}
	}
}

double DotProduct(double *a,double *b,int size){
	int i;
	double producto=0;
#pragma omp parallel for reduction(+:producto)
	for(i=0;i<size;i++){
		producto+=a[i]*b[i];
	}
	return producto;
}

double *DotProductVector(double e,double *v,int size){
	double *v_tem;
	int i;
	v_tem=(double*)malloc(sizeof(double)*size);
#pragma omp parallel for 
	for(i=0;i<size;i++){
		v_tem[i]=v[i]*e;
	}
	return v_tem;
}

void DotProductVector_E(double e,double *v,int size){
	int i;
#pragma omp parallel for 
	for(i=0;i<size;i++){
		v[i]=v[i]*e;
	}
}

double *SumVectors(double *a,double *b,int size){
	int i;
	double *suma=(double*)malloc(sizeof(double)*size);
#pragma omp parallel for 
	for(i=0;i<size;i++){
		suma[i]=a[i]+b[i];
	}
	return suma;
}


void CopyColumn(double *a,double *b,int size){
	int i;
#pragma omp parallel for 
	for(i=0;i<size;i++){
		a[i]=b[i];
	}
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

double** P_Process(double **K,double miu, double **M,int size){
	double **res= NULL;
	res = AllocateMatrixSpace(size);
	int i,j;
	#pragma omp parallel for private(j)
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			res[i][j]=((K[i][j])-(miu*M[i][j]));
		}
	}
	return res;
}

int LU_Method(double **A,int size){
	int k,sigma=0;
	double *U = AllocateVectorSpace(size);
	double *L= AllocateVectorSpace(size);
	double *a = AllocateVectorSpace(size);
	double *b = AllocateVectorSpace(size);
	double *c = AllocateVectorSpace(size);
	for(k=1;k<size;k++){
		a[k]= A[k-1][k];
	}
	for(k=0;k<size;k++){
		b[k]= A[k][k];
	}
	for(k=0;k<size-1;k++){
		c[k]= A[k+1][k];
	}
	U[0]=A[0][0];
	if(U[0]<0)
		sigma++;
	for (k=0;k<size-1;k++){
		L[k+1]=a[k+1]/U[k];
		U[k+1]=b[k+1]-(L[k+1]*c[k]);
		if(U[k+1]<0)
			sigma++;
	}
	return sigma;
}
double * AllocateVectorSpace(int size){
	double *m = (double*)malloc(sizeof(double)*size);
	return m;
}

double fuvx(double u,double x,double lambda){
	double ans=0.0;
		ans= (((-1.0)/(x*x*x))-lambda)*u;
	return ans;
}
