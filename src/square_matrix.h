#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#define EQUAL_VALUES 1
#define NOT_EQUAL_VALUES 0

void applyTranspose(double**, int);
void printSqrMat(const double**, int);
void printDelimSqrMat(const double**, char, int);
int areEqualsSqrMat(const double**, const double**, int);
double** getZeroSqrMat(int);
double** fillZeroSqrMat(double**, int);
double** freeSqrMat(double**, int);
double** getIdentitySqrMat(int);
double** fillIdentitySqrMat(double**, int);
double determinant(const double**, int);
double** comatrix(const double**, int);
double** applyComatrix(const double**, double**, int);
double** applyScalarMultiplication(double**, int, double);
double** getSqrMatInverse(const double**, int);
double** cloneSqrMat(const double**, int);
double** copySqrMat(const double**, double**, int);
double** getSqrMatMulti(const double**, const double**, int);
double** applySqrMatMulti(const double**, const double**, double**, int);

#endif
