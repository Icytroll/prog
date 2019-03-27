#include<stdio.h>
#include<stdlib.h>
#include"vector.h"
#include"matrix.h"

vector* vector_alloc(int n) {
	vector* b = (vector*)malloc(sizeof(vector));
	(*b).size = n;
	(*b).data = (double*)malloc(n*sizeof(double));
	return b;
}

void vector_set(vector* b, int i, double x){ 
	(*b).data[i] = x;
}

double vector_get(vector* b, int i) {
	return (*b).data[i];
}

vector* mv_mult(matrix* A, vector* b) {
	int n = A->size1, m = A->size2;
	
	vector* c = vector_alloc(n);
	
	double sum;
	for(int i=0;i<n;i++) {
		sum = 0;
		for(int j=0;j<m;j++)
			sum += matrix_get(A,i,j)*vector_get(b,j);
		vector_set(c,i,sum);
	}
	return c;
}

void vector_print(vector* b, char* s, FILE* stream) {
	fprintf(stream,"%s\n",s);
	int n = b->size;
	for(int i=0;i<n;i++) {
		fprintf(stream,"%7.3f \n",vector_get(b,i));
	}
	fprintf(stream,"\n");
}

void vector_free(vector* b) {
	free((*b).data);
	free(b);
}
