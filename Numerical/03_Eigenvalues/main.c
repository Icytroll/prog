#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void matrix_print(gsl_matrix* A, char* s, FILE* stream);

int jacobi_cyclic(gsl_matrix* V, gsl_matrix* D);
int jacobi_single_row(gsl_matrix* D, int p);

int main() {

/*--- A Jacobi diagonalization with cyclic sweeps ---*/
	
	int n = 4;
	double rnd;
	srand(time(NULL));
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);
	gsl_matrix* VD = gsl_matrix_calloc(n,n);
	gsl_matrix* VDVT = gsl_matrix_calloc(n,n);
	for(int i=0;i<n;i++) {
		for(int j=i;j<n;j++) {
			rnd = rand()/(double)RAND_MAX;
			gsl_matrix_set(A,i,j,rnd);
			gsl_matrix_set(A,j,i,rnd);
			gsl_matrix_set(D,i,j,rnd);
			gsl_matrix_set(D,j,i,rnd);
		}
	}
	
	printf("Initial A matrix ...\n");
	matrix_print(A,"A =",stdout);
	matrix_print(D,"D =",stdout);
	printf("Decomposing D ...\n");
	gsl_matrix_set_identity(V);
	jacobi_cyclic(V,D);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,D,0,VD);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,VD,V,0,VDVT);
	matrix_print(V,"V =",stdout);
	matrix_print(D,"D =",stdout);
	matrix_print(VDVT,"V*D*V' =",stdout);
	
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(D);
	gsl_matrix_free(VD);
	gsl_matrix_free(VDVT);

	clock_t t;
	double time_taken;
	FILE* A_time = fopen("cyclic_sweep_time.txt","w");
	for(int k=1;k<10;k++) {	
		n = pow(2,k);
		gsl_matrix* A = gsl_matrix_alloc(n,n);
		gsl_matrix* V = gsl_matrix_alloc(n,n);
		gsl_matrix* D = gsl_matrix_alloc(n,n);
		for(int i=0;i<n;i++) {
			for(int j=i;j<n;j++) {
				rnd = rand()/(double)RAND_MAX;
				gsl_matrix_set(A,i,j,rnd);
				gsl_matrix_set(A,j,i,rnd);
			}
		}
	
		gsl_matrix_memcpy(D,A);	
		gsl_matrix_set_identity(V);
		
		t = clock();
		jacobi_cyclic(V,D);
		t = clock() - t;
		time_taken = ((double)t)/CLOCKS_PER_SEC;
		fprintf(A_time,"%d %f\n",n,time_taken);
		
		gsl_matrix_free(A);
		gsl_matrix_free(V);
		gsl_matrix_free(D);
	}

	fclose(A_time);
	
/*--- B Jacobi diagonalization eigenvalue-by-eigenvalue ---*/
		
	n = 5;
	A = gsl_matrix_alloc(n,n);
	D = gsl_matrix_alloc(n,n);
	V = gsl_matrix_alloc(n,n);
	for(int i=0;i<n;i++) {
		for(int j=i;j<n;j++) {
			rnd = rand()/(double)RAND_MAX;
			gsl_matrix_set(A,i,j,rnd);
			gsl_matrix_set(A,j,i,rnd);
		}
	}
	gsl_matrix_memcpy(D,A);
	
	int eig = 3;
	printf("\n\nCalculating the first %d eigenvalues of D, one at a time ...\n",eig);
	matrix_print(D,"D = ",stdout);
	for(int p=0;p<eig;p++) {
		jacobi_single_row(D,p);
		matrix_print(D,"D =",stdout);
	}
	
	printf("\n\n");
	gsl_matrix_memcpy(D,A);
	int tot_rot;
	tot_rot = jacobi_cyclic(V,D);
	printf("Number of rotations in D for full cyclic sweep: %d\n",tot_rot);
	tot_rot = 0;
	gsl_matrix_memcpy(D,A);
	for(int p=0;p<n;p++)
		tot_rot += jacobi_single_row(D,p);
	printf("Total number of rotations in D for single row sweeps: %d\n",tot_rot);
	
/*--- C Classic Jacobi eigenvalue algorithm ---*/
	
	
}
