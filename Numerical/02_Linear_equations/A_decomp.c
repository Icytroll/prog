#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"

void qr_gs_decomp(matrix* Q, matrix* R);

void A_decomp(FILE * Astream) {
	fprintf(Astream,"Testing decomposing function on tall matrix ...\n\n");
	
	int n = 6, m = 4;
	matrix* A = matrix_alloc(n,m);
	matrix* Q = matrix_alloc(n,m);
	matrix* R = matrix_alloc(m,m);
	
	// Fill A and Q with the same random numbers between 0 and 1.
	// A stays as it is, while Q will be changed in qr_gs_decomp.
	double rnd;
	srand(time(NULL));
	for(int i=0;i<n;i++) {
		for(int j=0;j<m;j++) {
			rnd = rand()/(double)RAND_MAX;	
			matrix_set(A,i,j,rnd);
			matrix_set(Q,i,j,rnd);
		}
	}
	matrix_print(A,"A =",Astream);
	fprintf(Astream,"Making Q a copy of A:\n");
	matrix_print(Q,"Q =",Astream);
	qr_gs_decomp(Q,R);

	// Check that everything is as it should be
	fprintf(Astream,"After decomposition:\n");
	matrix_print(Q,"Q =",Astream);
	matrix_print(R,"R =",Astream);
	
	matrix* QT = matrix_transpose(Q);
	matrix* QT_Q = matrix_alloc(n,m);
	matrix_mult(QT,Q,QT_Q);
	
	matrix_print(QT_Q,"QT*Q =",Astream);
	matrix* QR = matrix_alloc(n,m);
	matrix_mult(Q,R,QR);
	matrix_print(QR,"Q*R =",Astream);
	
	// Free memory
	matrix_free(A);
	matrix_free(Q);
	matrix_free(R);
	matrix_free(QT);
	matrix_free(QT_Q);
	matrix_free(QR);

}


