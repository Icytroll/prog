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
