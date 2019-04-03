#include<stdio.h>
#include<math.h>

#ifndef HAVE_VECTOR_H /* this is necessary for multiple includes */
typedef struct {int size; double *data;} vector;

vector* vector_alloc(int n);
void vector_set(vector* b, int i, double x);
double vector_get(vector* b, int i);
double norm(vector* b);
void vector_mult_all(vector* b, double x);
void vector_add(vector* a, vector* b, vector* c);
void vector_sub(vector* a, vector* b, vector* c);
void vector_print(vector* b, char* s, FILE *stream);
void vector_free(vector* b);

#define HAVE_vector_H
#endif
