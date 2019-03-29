#include"matrix.h"
#include"vector.h"

void mv_mult(matrix* A, vector* b, vector* c);

void bi_inverse(matrix* UT, matrix* B, matrix* V, matrix* invA) {
    // A*x_i = e_i --> U*B*V'*x_i = e_i --> B*V'*x_i = UTe_i --> B*V'*x_i = UT_i
    // With e_i being the ith coloumn of an identity matrix, x_i
    // will be the ith coloumn of the inverse of A.

    int n = invA->size1;
	vector* y = vector_alloc(n);
    vector* x = vector_alloc(n);
	double s;
    for(int i=0;i<n;i++) {
        for(int j=n-1;j>=0;j--) {
            s = matrix_get(UT,j,i);
            for(int k=j+1;k<n;k++) {
                s -= matrix_get(B,j,k)*vector_get(y,k);
            }
            vector_set(y,j,s/matrix_get(B,j,j));
        }
		mv_mult(V,y,x);
		for(int j=0;j<n;j++) matrix_set(invA,j,i,vector_get(x,j));
    }
	vector_free(y);
	vector_free(x);
}
