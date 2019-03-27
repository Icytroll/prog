#include<math.h>
#include"matrix.h"

void bi_decomp(matrix* A, matrix* U, matrix* B, matrix* V) {
	
	int n = A->size1;
	double a,b = 0,sum;
	
	// Set first column of V as a unit vector
	matrix_set(V,0,0,1);
	

	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++){
			sum = 0;
			for(int k=0;k<n;k++) sum += matrix_get(A,j,k)*matrix_get(V,k,i);
			if(i!=0) sum -= b*matrix_get(U,j,i);
			matrix_set(U,j,i,sum);
		}
		sum = 0;
		for(int j=0;j<n;j++) sum += matrix_get(U,j,i)*matrix_get(U,j,i);
		a = sqrt(sum);
		matrix_set(B,i,i,a);
		for(int j=0;j<n;j++) matrix_set(U,j,i,matrix_get(U,j,i)/a);
		if(i == n) break;

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
