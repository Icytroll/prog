#include"matrix.h"

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_inverse(matrix* Q, matrix* R, matrix* invA);

void inverse(matrix* A, matrix* invA) {
	
	int n = A->size1;
	matrix* Q = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			matrix_set(Q,i,j,matrix_get(A,i,j));

	qr_gs_decomp(Q,R);
	qr_gs_inverse(Q,R,invA);
}
