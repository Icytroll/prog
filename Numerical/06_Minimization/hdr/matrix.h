#include<stdio.h>
#include<math.h>

#ifndef HAVE_MATRIX_H /* this is necessary for multiple includes */
typedef struct {int size1, size2; double *data;} matrix;

matrix* matrix_alloc(int n, int m);
void matrix_set(matrix* A, int i, int j, double x);
double matrix_get(matrix* A, int i, int j);
void matrix_add(matrix* A, matrix* B, matrix* C);
void matrix_transpose(matrix* A, matrix* AT);
void matrix_set_identity(matrix* A);
void matrix_mult(matrix* A, matrix* B, matrix* C);
void matrix_print(matrix* A, char* s, FILE* stream);
void matrix_free(matrix* A);

#define HAVE_MATRIX_H
#endif
