#ifndef QRFACT_H
#include <stdio.h>
#define QRFACT_H


double** AllocateMatrixSpace(int);
double** InitializingA(double, double, double);
void Transpose(double**);
double** StartM(int);
double** StartQ(double, double, int);
void QR_Method(double**);
double** Multiplication(double**, double**);
void CopyMatrix(double**, double**);
double DotProduct(double*,double*);
void DotProductVector_E(double,double*);
double* DotProductVector(double,double*);
double* SumVectors(double*,double*);
void PrintMatrix(double**,int);
void CopyColumn(double*,double*);
void WriteToFile(FILE*, double**,int);


#endif
