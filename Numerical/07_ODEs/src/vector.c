#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"vector.h"
#include"matrix.h"

vector* vector_alloc(int n) {
	vector* b = (vector*)malloc(sizeof(vector));
	(*b).size = n;
	(*b).data = (double*)calloc(n,sizeof(double));
	return b;
}

void vector_set(vector* b, int i, double x){ 
	(*b).data[i] = x;
}

void vector_set_all(vector* b, double x){
	int n = b->size;
	for(int i=0;i<n;i++) vector_set(b,i,x);
}

double vector_get(vector* b, int i) {
	return (*b).data[i];
}

double norm(vector* b) {
	int n = b->size;
	double sum = 0;
	for(int i=0;i<n;i++) sum += vector_get(b,i)*vector_get(b,i);
	return sqrt(sum);
}

double vector_inner(vector* a, vector* b) {
	int n = a->size;
	double sum = 0;
	for(int i=0;i<n;i++) sum += vector_get(a,i)*vector_get(b,i);
	return sum;
}

void vector_mult_all(vector* b, double x) {
	int n = b->size;
	for(int i=0;i<n;i++) vector_set(b,i,vector_get(b,i)*x);
}

void vector_add(vector* a, vector* b, vector* c) {
	int n = a->size;
	for(int i=0;i<n;i++)
		vector_set(c,i,vector_get(a,i)+vector_get(b,i));
}

void vector_sub(vector* a, vector* b, vector* c) {
	int n = a->size;
	for(int i=0;i<n;i++)
		vector_set(c,i,vector_get(a,i)-vector_get(b,i));
}

void vector_print(vector* b, char* s, FILE* stream) {
	fprintf(stream,"%s\n",s);
	int n = b->size;
	for(int i=0;i<n;i++) {
		fprintf(stream,"%14.10f \n",vector_get(b,i));
	}
	fprintf(stream,"\n");
}

void vector_free(vector* b) {
	free((*b).data);
	free(b);
}
