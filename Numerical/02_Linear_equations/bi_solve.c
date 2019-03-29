#include"matrix.h"
#include"vector.h"

void mv_mult(matrix* A, vector* b, vector* c);

void bi_solve(matrix* UT, matrix* B, matrix* V, vector* x, vector* b) {
	int n = B->size1;
	vector* UTb = vector_alloc(n);
	mv_mult(UT,b,UTb);
	vector* y = vector_alloc(n);
	double s;
	for(int i=n-1;i>=0;i--) {
		s = vector_get(UTb,i);
		for(int j=i+1;j<n;j++) s -= matrix_get(B,i,j)*vector_get(y,j);
		vector_set(y,i,s/matrix_get(B,i,i));
	}
	mv_mult(V,y,x);
	vector_free(UTb);
	vector_free(y);
}
