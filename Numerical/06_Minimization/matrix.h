#include<stdio.h>

#ifndef HAVE_MATRIX_H /* this is necessary for multiple includes */
typedef struct {int size1, size2; double *data;} matrix;

matrix* matrix_alloc(int n, int m);
void matrix_set(matrix* A, int i, int j, double x);
double matrix_get(matrix* A, int i, int j);
matrix* matrix_transpose(matrix* A);
matrix* matrix_mult(matrix* A, matrix* B);
void matrix_print(matrix* A, char* s, FILE* stream);
void matrix_free(matrix* A);

#define HAVE_MATRIX_H
#endif
