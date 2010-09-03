#include "qrfact.h"
#include "square_matrix.h"
#include "magnetic_monopole.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

void test_mm_numerov(int size) {
  // Numerov
  double start=0.001;
  double final=3.1416;
  int    qr_iterations = 100;

  printf("Magnetic Monopole Spectrum (Numerov Method)\n");
  printf("Start:\t%lf\n",start);
  printf("Final:\t%lf\n",final);
  printf("Size:\t%d\n",size);
  printf("CPUs:\t%d\n",omp_get_num_procs());
  printf("Iters:\t%d\n",qr_iterations);

  double** mm_num = mm_numerov(start, final, qr_iterations, size);
  printf("Matrix A\n");
  printSqrMat(mm_num,size);
  freeSqrMat(mm_num,size);
}

void test_mm_rungekutta(int size) {
  // Runge-Kutta
  double delta  = 0.001;
  double start  = 0.001;
  double final  = 3.1416;
  double h      = (final-start)/size;

  printf("Magnetic Monopole Spectrum (Runge-Kutta Method)\n");
  printf("Start:\t%lf\n",start);
  printf("Final:\t%lf\n",final);
  printf("Size:\t%d\n",size);
  printf("CPUs:\t%d\n",omp_get_num_procs());
  printf("H:\t%lf\n",h);

  mm_rungekutta(start,final,size,delta);
}

void test_mm_lambda(int size) {
  // Lambda
  double delta  = 0.001;
  double start  = 0.001;
  double final  = 3.1416;
  double h      = (final-start)/size;
  double miu   = -15.0;

  printf("Magnetic Monopole Spectrum (Lambda Method)\n");
  printf("Start:\t%lf\n",start);
  printf("Final:\t%lf\n",final);
  printf("Size:\t%d\n",size);
  printf("CPUs:\t%d\n",omp_get_num_procs());
  printf("H:\t%lf\n",h);
  printf("Delta:\t%lf\n",delta);
  printf("Miu:\t%lf\n",miu);

  mm_lambda(start,final,size,delta,miu);
}

int main(int argc, char** argv) {

  if(argc!=2) {
    fprintf(stderr, "Usage: %s <size>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  omp_set_num_threads(omp_get_num_procs()*2);

  int size = strtol(argv[1], (char **)NULL, 10);

  test_mm_numerov(size);
  test_mm_lambda(size);
  test_mm_rungekutta(size);
 
  return 0;
}
