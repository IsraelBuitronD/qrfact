#ifndef QRFACT_H
#include <stdio.h>
#define QRFACT_H


double** AllocateMatrixSpace();
double** InitializingA(double, double, double);
void Transpose(double**);
double **StartM();
double **StartQ(double,double, double);
double** K_Process(double **);
void QR_Method(double**);
double** Multiplication(double**, double**);
void CopyMatrix(double**, double**);
double DotProduct(double*,double*);
void DotProductVector_E(double,double*);
double* DotProductVector(double,double*);
double* SumVectors(double*,double*);
void PrintMatrix(double**);
void CopyColumn(double*,double*);
void WriteToFile(FILE*, double**,int);


#endif
