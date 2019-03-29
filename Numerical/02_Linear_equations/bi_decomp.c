#include<math.h>
#include<assert.h>
#include"matrix.h"

void bi_decomp(matrix* A, matrix* U, matrix* B, matrix* V) {
	assert(A->size1==U->size1);
	assert(A->size1==B->size1);
	assert(A->size1==V->size1);
	assert(A->size2==A->size1);
	assert(U->size2==A->size1);
	assert(B->size2==A->size1);
	assert(V->size2==A->size1);
	
	int n = A->size1;
	double a,b = 0,sum;
	
	// Set first column of V as a unit vector
	matrix_set(V,0,0,1);
	

	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++){
			sum = 0;
			for(int k=0;k<n;k++) sum += matrix_get(A,j,k)*matrix_get(V,k,i);
			if(i!=0) sum -= b*matrix_get(U,j,i-1);
			matrix_set(U,j,i,sum);
		}
		sum = 0;
		for(int j=0;j<n;j++) sum += matrix_get(U,j,i)*matrix_get(U,j,i);
		a = sqrt(sum);
		matrix_set(B,i,i,a);
		for(int j=0;j<n;j++) matrix_set(U,j,i,matrix_get(U,j,i)/a);
		if(i == n-1) break;

		for(int j=0;j<n;j++) {
			sum = 0;
			for(int k=0;k<n;k++) sum+= matrix_get(A,k,j)*matrix_get(U,k,i);
			matrix_set(V,j,i+1,sum - a*matrix_get(V,j,i));
		}
		sum = 0;
		for(int j=0;j<n;j++) sum += matrix_get(V,j,i+1)*matrix_get(V,j,i+1);
		b = sqrt(sum);
		matrix_set(B,i,i+1,b);
		for(int j=0;j<n;j++) matrix_set(V,j,i+1,matrix_get(V,j,i+1)/b);
	}
}
