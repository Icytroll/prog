#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"
#include"vector.h"

void mv_mult(matrix* A, vector* b, vector* c);
void bi_decomp(matrix* A, matrix* U, matrix* B, matrix* V);
double bi_det(matrix* B);
void bi_solve(matrix* UT, matrix* B, matrix* V, vector* x, vector* b);
void bi_inverse(matrix* UT, matrix* B, matrix* V, matrix* invA);

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

	int n = 6, m = 6;
	matrix* A = matrix_alloc(n,m);
	matrix* U = matrix_alloc(n,m);
	matrix* B = matrix_alloc(n,m);
	matrix* V = matrix_alloc(n,m);
	vector* b = vector_alloc(n);
	vector* x = vector_alloc(n);

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
	
	// Decomposing A
	bi_decomp(A,U,B,V);
	matrix_print(U,"U =",Cstream);
	matrix_print(B,"B =",Cstream);
	matrix_print(V,"V =",Cstream);
	matrix* VT = matrix_transpose(V);
	matrix* UT = matrix_transpose(U);
	matrix* UB = matrix_alloc(n,n);
	matrix_mult(U,B,UB);
	matrix* UBVT = matrix_alloc(n,n);
	matrix_mult(UB,VT,UBVT);
	matrix_print(UBVT,"U*B*V' =",Cstream);
	
	// Solving linear system of equations
	fprintf(Cstream,"\nSolving U*B*V'*x = b ...\n");
	bi_solve(UT,B,V,x,b);
	vector_print(x,"x =",Cstream);
	vector* UBVTx = vector_alloc(n);
	vector* Ax = vector_alloc(n);
	mv_mult(A,x,Ax);
	mv_mult(UBVT,x,UBVTx);
	vector_print(Ax,"A*x =",Cstream);
	vector_print(UBVTx,"U*B*V'*x =",Cstream);
	
	// Calculate determinant
	fprintf(Cstream,"\nDeterminant of the matrix B: %g\n",-bi_det(B));
	

	// Calculate inverse of A
	fprintf(Cstream,"\n\nInverting A via U, B and V ...\n");
	matrix* invA = matrix_alloc(n,m);
	bi_inverse(UT,B,V,invA);
	matrix_print(invA,"A^-1 =",Cstream);
	matrix* AinvA = matrix_alloc(n,n);
	matrix_mult(A,invA,AinvA);
	matrix_print(AinvA,"A*invA =",Cstream);
	
	// Clean up
	matrix_free(A);
	matrix_free(U);
	matrix_free(B);
	matrix_free(V);
	vector_free(x);
	vector_free(b);
	matrix_free(VT);
	matrix_free(UT);
	matrix_free(UB);
	matrix_free(UBVT);
	vector_free(UBVTx);
	vector_free(Ax);
	matrix_free(invA);
	matrix_free(AinvA);
}
