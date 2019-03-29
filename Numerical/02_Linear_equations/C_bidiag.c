#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"
#include"vector.h"

void bi_decomp(matrix* A, matrix* U, matrix* B, matrix* V);


void C_bidiag(FILE * Cstream) {
	
	/* Determinant
	double detB = 1;
	for(int i=0;i<n;i++)
		detB *= matrix_get(B,i,i);
	fprintf(Cstream,"Determinant of A = %g\n\n",detB);
	
	Ax = b
	UBV'x = b
	BV'x = U'b
	V'x = y
	By = U'b --> y
	V'x = y --> x
	
	Ax = ei
	UBV'x = ei
	BV'x = U'ei
	By = U'b --> y
	V'x = y --> x
	
	*/

	fprintf(Cstream,"Testing bidiagonalization ...\n\n");


	printf("bidiag called...\n");
	
	int n = 4, m = 4;
	matrix* A = matrix_alloc(n,m);
	matrix* U = matrix_alloc(n,m);
	matrix* B = matrix_alloc(n,m);
	matrix* V = matrix_alloc(n,m);
	vector* b = vector_alloc(n);
	vector* x = vector_alloc(n);

	printf("bidiag: all allocated ...\n");
	
	// Initialize U,B and V to zero
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			matrix_set(U,i,j,0);
			matrix_set(B,i,j,0);
			matrix_set(V,i,j,0);
		}
	}

	double rnd;
	srand(time(NULL));
	for(int i=0;i<n;i++) {
		for(int j=i;j<m;j++) {
			// Fill A symmetrically with random numbers from 0 to 1
			rnd = rand()/(double)RAND_MAX;	
			matrix_set(A,i,j,rnd);
			matrix_set(A,j,i,rnd);
		}
		// Fill b with random numbers from 0 to 1
		rnd = rand()/(double)RAND_MAX;	
		vector_set(b,i,rnd);
	}
	
	matrix_print(A,"A =",Cstream);
	vector_print(b,"b =",Cstream);

	printf("bidiag: calling bi_decomp...\n");
	bi_decomp(A,U,B,V);
	printf("bidiag: bi_decomp exited...\n");
	matrix_print(U,"U =",Cstream);
	matrix_print(B,"B =",Cstream);
	matrix_print(V,"V =",Cstream);
	
	matrix* UT = matrix_transpose(U);
	printf("bidiag: calling transpose...\n");
	matrix* VT = matrix_transpose(V);
	printf("This is ok\n");
	matrix* UTU = matrix_mult(UT,U);
	matrix* VTV = matrix_mult(VT,V);
	matrix_print(UT,"U' =",Cstream);
	matrix_print(VT,"V' =",Cstream);
	matrix_print(UTU,"U'*U =",Cstream);
	matrix_print(VTV,"V'*V =",Cstream);
	
	printf("calling matrix_mult...\n");
	matrix* UB = matrix_mult(U,B);
	printf("after matrix_mult: This is not ok???\n");
	
	
	matrix* UBVT = matrix_mult(UB,VT);
	matrix_print(UBVT,"U*B*V' =",Cstream);
	
	/*
	bi_solve(U,B,V,x,b);
	double detA = bi_det(B);

	matrix* invA = matrix_alloc(n,m);

	bi_inverse(U,B,V,invA);

	*/
}
