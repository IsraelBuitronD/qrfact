#include "qrfact.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

int SIZE;
double **q,**r;

int main(int argc, const char * argv[]){
	FILE *M_A,*M_Q,*M_R;
	M_A =fopen("Matriz_AQ.txt", "w"); 
	M_Q =fopen("Matriz_Q.txt", "w");
	M_R= fopen("Matriz_R.txt","w");
	SIZE=0;
	int i;
	double start=0.0, final=0.0;
	double **A;
	printf("QR iterative process\n");
	//Getting Configuration Values
	start=(double)  atoi(argv[1]);
	final=(double)  atoi(argv[2]);
	SIZE= atoi(argv[3]);
	printf("Start Point: %lf\n",start);
    printf("Final Point:  %lf\n ", final);
	printf("Partition Size: %d\n", SIZE);
	printf("CPUs Number: %d\n",omp_get_num_procs());
	omp_set_num_threads(omp_get_num_procs()*2);
	q=AllocateMatrixSpace();
	r=AllocateMatrixSpace();
    A=InitializingA(start,final,SIZE);
	// Iterative proccess 
	for (i=0;i<80;i++){
		QR_Method(A);
		//Transpose(q);
		A=Multiplication(q,A);
		Transpose(q);
		A=Multiplication(q,A);
	}
	//WriteToFile(M_A,A,SIZE);
	WriteToFile(M_Q,q,SIZE);
	WriteToFile(M_R,r,SIZE);
	
}

double **AllocateMatrixSpace(){
	int i;
	double **m;
	m=(double**)malloc(sizeof(double)*SIZE);
	for(i=0;i<SIZE;i++){
		m[i]=(double*)malloc(sizeof(double)*SIZE);
	}
	return m;
}

double **InitializingA(double StartPoint, double FinalPoint, double PartitionNumber){
	double **result;
	result=AllocateMatrixSpace();
	int i,j;
	double PartitionSize=(FinalPoint-StartPoint)/PartitionNumber;
	double T,Q,div,actual=StartPoint;
	div=1/(PartitionSize*PartitionSize);
	for(i=0;i<PartitionNumber;i++){
		for(j=0;j<PartitionNumber;j++){
			if(i==j){
				T=2;
				actual=StartPoint+PartitionSize*i;
				Q=(actual-0.5)*(actual-0.5); // Function Implementation (X-0.5)
			}
			else{
				if( (i<PartitionNumber-1&&(i==j+1||i==j-1))||(i>0&&(i==j-1||i==j+1 ) )){
					T=-1;
				}
				else
					T=0;
				Q=0;
			}
			result[i][j]=(div*T)+Q;
		}
	}
	return result;
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

void QR_Method(double **a){
	int i,j;
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


double **Multiplication(double **matriz1,double **matriz2){
	int i,j,k;
	double **matrizR=AllocateMatrixSpace();
	for(i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			for(k=0;k<SIZE;k++){
				matrizR[i][j]+=matriz1[i][k]*matriz2[k][j];
			}
		}
	}
	return matrizR;
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

