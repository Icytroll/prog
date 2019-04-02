#include<math.h>
#include"matrix.h"
#include"vector.h"

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);
void inverse(matrix* A, matrix* invA);

void lsfit(
	int n, double funs(int,double),
	double* x, double* y, double* dy,
	vector* c, matrix* S)
{
	int m = c->size;
	matrix* A = matrix_alloc(n,m);
    vector* b = vector_alloc(n);
    matrix* Q = matrix_alloc(n,m);
    matrix* R = matrix_alloc(m,m);

    // Build matrix A (and Q) and vector b from data points
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            matrix_set(A,i,j,funs(j,x[i])/dy[i]);
            matrix_set(Q,i,j,matrix_get(A,i,j));
        }
        vector_set(b,i,y[i]/dy[i]);
    }
	

	// Decompse A (Q) into Q and R
	qr_gs_decomp(Q,R);

	// Calculate coefficients c
	qr_gs_solve(Q,R,b,c);
    
	// Calculate Covariance matrix
	matrix* RT = matrix_alloc(m,m);
    matrix_transpose(R,RT);
    matrix* RTR = matrix_alloc(m,m);
    matrix_mult(RT,R,RTR);
    inverse(RTR,S);
	
	//matrix_print(A,"A =",stdout);
	matrix_free(A);
	matrix_free(Q);
	matrix_free(R);
	vector_free(b);
	matrix_free(RT);
	matrix_free(RTR);	
}
