#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

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
double** applyScalarMultiplication(double**, int, double);
double** getSqrMatInverse(double**, int);
double** cloneSqrMat(double **, int);

#endif
