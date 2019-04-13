#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"

matrix* matrix_alloc(int n, int m) {
	matrix* A=(matrix*)malloc(sizeof(matrix));
	(*A).size1=n; (*A).size2=m;
	(*A).data=(double*)calloc(n*m,sizeof(double));
	return A;
}

void matrix_set(matrix* A, int i, int j, double x){ 
	 (*A).data[i+j*(*A).size1] = x;
}

double matrix_get(matrix* A, int i, int j) {
	return (*A).data[i+j*(*A).size1];
}

void matrix_add(matrix* A, matrix* B, matrix* C) {
	int n = A->size1, m = A->size2;
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			matrix_set(C,i,j,matrix_get(A,i,j)+matrix_get(B,i,j));
}

void matrix_transpose(matrix* A, matrix* AT) {
	int n = A->size2, m = A->size1;
	for(int i=0;i<n;i++) {
		for(int j=0;j<m;j++) {
			matrix_set(AT,i,j,matrix_get(A,j,i));
		}
	}
}

void matrix_set_identity(matrix* A) {
	int n = A->size1;
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			if (i==j) matrix_set(A,i,j,1);
			else matrix_set(A,i,j,0);
		}
	}
}

#include<assert.h>
void matrix_mult(matrix* A, matrix* B, matrix* C) {
	assert(A->size2==B->size1);
	int n = A->size1, m = B->size2;
	int l = A->size2;
	double sum;
	for(int i=0;i<n;i++) {
		for(int j=0;j<m;j++) {
			sum = 0;
			for(int k=0;k<l;k++)
				sum += matrix_get(A,i,k)*matrix_get(B,k,j);
			matrix_set(C,i,j,sum);
		}
	}
}

void matrix_print(matrix* A, char* s, FILE* stream) {
	fprintf(stream,"%s\n",s);
	int n = A->size1;
	int m = A->size2;
	for(int i=0;i<n;i++) {
		for(int j=0;j<m;j++) {
			fprintf(stream,"%7.3f ",matrix_get(A,i,j));
		}
		fprintf(stream,"\n");
	}
	fprintf(stream,"\n");
}

void matrix_free(matrix* A) {
	free(A->data);
	free(A);
}
