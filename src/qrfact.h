#ifndef QRFACT_H
#include <stdio.h>
#define QRFACT_H


double** AllocateMatrixSpace(int);
double** InitializingA(double, double, double);
double** StartM(double***, int);
double** StartQ(double, double, double*, int);
double** K_Process(double**, double**, double**, double*, int);
double** QR_Method(double**,int);
double** Multiplication(double**, double**, int);
void CopyMatrix(double**, double**,int);
double DotProduct(double*,double*,int);
void DotProductVector_E(double,double*,int);
double* DotProductVector(double,double*,int);
double* SumVectors(double*,double*,int);
void PrintMatrix(double**,int);
void CopyColumn(double*,double*,int);
void WriteToFile(FILE*, double**,int);
double **P_Process(double**,double,double**,int size);
int LU_Method(double **,int);
double * AllocateVectorSpace(int);

#endif
