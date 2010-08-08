#ifndef QRFACT_H
#define QRFACT_H

void qr(double**);
double** pedirMemoriaMatriz();
double* pedirMemoriaVector();
void inicializaA(double**);
void copiarMatriz(double**, double**);
void transpuesta(double**);
void imprime(double**);
double productoPunto(double*,double*);
void productoEscalarVector1(double,double*);
double* productoEscalarVector(double,double*);
double* sumaVectores(double*,double*);
void copiarColumna(double*,double*);
double** crearA(double, double, double);
double** multiplica(double**, double**);

#define SIZE 10

#endif
