#ifndef MATRIX_H
#define MATRIX_H

#define EQUAL_VALUES 1
#define NOT_EQUAL_VALUES 0

void applyTranspose(double**, int);
//double** getTranspose(double**, int);
void printSqrMat(double**, int);
void printDelimSqrMat(double**, char, int);
int areEqualsSqrMat(double**, double**, int);
double** getZeroSqrMat(int);
double** freeSqrMat(double**, int);
double** getIdentitySqrMat(int);
double determinant(double**, int);
double** comatrix(double **, int);

#endif
