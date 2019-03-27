#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"
#include"vector.h"

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* QT, matrix* R, vector* b, vector* x);

void A_solve(FILE * Astream) {
	
	fprintf(Astream,"Solving linear system of equations ...\n\n");
	int n = 6;
	matrix* A = matrix_alloc(n,n);
	matrix* Q = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	vector* b = vector_alloc(n);
	vector* x = vector_alloc(n);
	
	// Fill A and Q with the same random numbers between 0 and 1.
	// A stays as it is, while Q will be changed in qr_gs_decomp.
	double rnd;
	srand(time(NULL));
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			rnd = rand()/(double)RAND_MAX;	
			matrix_set(A,i,j,rnd);
			matrix_set(Q,i,j,rnd);
		}
	}

	// Fill out b
	for(int i=0;i<n;i++)
		vector_set(b,i,rand()/(double)RAND_MAX);
	
	// Factorize A into Q and R, and solve Rx=QTb for x
	qr_gs_decomp(Q,R);
	matrix* QT = matrix_transpose(Q);
	qr_gs_solve(QT,R,b,x);

	// Check solution
	vector* Ax = vector_alloc(n);
	double sum;
	for(int i=0;i<n;i++) {
		sum = 0;
		for(int j=0;j<n;j++) sum += matrix_get(A,i,j)*vector_get(x,j);
		vector_set(Ax,i,sum);
	}
	
	// Print out everything
	matrix_print(A,"A =",Astream);
	vector_print(b,"b =",Astream);
	matrix_print(Q,"Q =",Astream);
	matrix_print(R,"R =",Astream);
	vector_print(x,"x =",Astream);
	vector_print(Ax,"A*x =",Astream);
	
	// Free memory
	matrix_free(A);
	vector_free(b);
	matrix_free(Q);
	matrix_free(R);
	vector_free(x);
	matrix_free(QT);
	vector_free(Ax);


}
