#ifndef MAGNETIC_MONOPOLE_H
#define MAGNETIC_MONOPOLE_H

double** mm_numerov(double, double, int, int);
double mm_rungekutta(double, double, int, double);
double mm_lambda(double, double, int, double, double);
double function_uvx(double, double, double);
double** startQ(double, double, int);
double** kProcess(double**, double**, double**, double, int);
double** qrMethod(double**, int);
int luMethod(double**, int);
double** pProcess(double**, double, double**, int);
double** fillT(int);
double** startM(double**, int);

#endif
