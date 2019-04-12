#include"matrix.h"
#include"vector.h"

void mv_mult(int transposeA, matrix* A, vector* b, vector* c);

void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x) {
	
	// Apply QT to b, send the output to x
	mv_mult(1,Q,b,x);
	
	// Perform in-place back-substitution on x
	int n = R->size1;
	for(int i=n-1;i>=0;i--) {
		double s = vector_get(x,i);
		for(int j=i+1;j<n;j++) s -= matrix_get(R,i,j)*vector_get(x,j);
		vector_set(x,i,s/matrix_get(R,i,i));
	}
}
