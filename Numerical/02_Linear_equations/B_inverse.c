#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_inverse(matrix* Q, matrix* R, matrix* B);

void B_inverse(FILE * Bstream) {
	
	fprintf(Bstream,"Inverting square matrix A ...\n\n");
	int n = 6;
	matrix* A = matrix_alloc(n,n);
	matrix* B = matrix_alloc(n,n);
	matrix* Q = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	
	// Initilize R to zeros
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			matrix_set(R,i,j,0);

	// Fill A and Q with the same random integers between 0 and 9.
	// A stays as it is, while Q will be changed in qr_gs_decomp.
	int rnd;
	srand(time(NULL));
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			rnd = round(rand()/(double)RAND_MAX*9);	
			matrix_set(A,i,j,rnd);
			matrix_set(Q,i,j,rnd);
		}
	}

	// Factorize A into Q and R, and use them to invert the matrix A into B
	qr_gs_decomp(Q,R);
	qr_gs_inverse(Q,R,B);

	matrix* AB = matrix_mult(A,B);
	matrix* QR = matrix_mult(Q,R);

	matrix_print(A,"A =",Bstream);
	matrix_print(Q,"Q =",Bstream);
	matrix_print(R,"R =",Bstream);
	matrix_print(QR,"Q*R =",Bstream);
	matrix_print(B,"B =",Bstream);
	matrix_print(AB,"A*B =",Bstream);

	matrix_free(A);
	matrix_free(Q);
	matrix_free(R);
	matrix_free(B);
	matrix_free(QR);
	matrix_free(AB);

}
