#include<stdio.h>
#include<stdlib.h>

void A_decomp(FILE * Astream);
void A_solve(FILE * Astream);
void B_inverse(FILE * Bstream);
void C_bidiag(FILE * Cstream);

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
	
	// Invert a random matrix by solving sets of linear equations
	B_inverse(Bstream);
	
	fclose(Bstream);

/*----- Golub-Kahan-Lanczos bidiagonalization -----*/
	
	FILE * Cstream = fopen("C.txt","w");
	
	C_bidiag(Cstream);	
	
	fclose(Cstream);
	
	return 0;	
}
