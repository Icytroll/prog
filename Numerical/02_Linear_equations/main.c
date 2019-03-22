#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"
#include"vector.h"

void A_decomp(FILE * Astream);
void A_solve(FILE * Astream);

//void qr_gs_inverse (const matrix* Q, const matrix* R, matrix* B);


int main() {
	
/*----- QR-decomposition by modified Gram-Schmidt orthogonalization -----*/
	
	FILE * Astream = fopen("A.txt","w");
	
	// Decompose a random tall matrix
	A_decomp(Astream);
	
	// Solve a random system of linear equations
	A_solve(Astream);
	
	fclose(Astream);
	
/*----- Matrix inverse by Gram-Schmidt QR factorization -----*/
	
	FILE * Bstream = fopen("B.txt","w");
	
	fclose(Bstream);

/*----- Golub-Kahan-Lanczos bidiagonalization -----*/
	
	FILE * Cstream fopen("C.txt","w");
	
	
	
	fclose(Cstream);
	
	return 0;	
}
