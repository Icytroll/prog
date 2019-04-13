#include"vector.h"
#include"matrix.h"

void vector_outer(vector* a, vector* b, matrix* C) {
	int n = a->size;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			matrix_set(C,i,j,vector_get(a,i)*vector_get(b,j));
}
