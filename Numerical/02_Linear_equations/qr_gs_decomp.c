#include<math.h>
#include"matrix.h"

void qr_gs_decomp(matrix* A, matrix* R) {
	int n = A->size1, m = A->size2;
    double R_ii, R_ij, a_jk, a_ik;
	for(int i=0;i<m;i++) {
		
		R_ii = 0;
        for(int k=0;k<n;k++) {
			a_ik = matrix_get(A,k,i);
            R_ii += a_ik*a_ik;
		}
        R_ii = sqrt(R_ii);
		matrix_set(R,i,i,R_ii);
		for(int k=0;k<n;k++)
			matrix_set(A,k,i,matrix_get(A,k,i)/R_ii);
		for(int j=i+1;j<m;j++) {
			
			R_ij = 0;
			for(int k=0;k<n;k++)
				R_ij += matrix_get(A,k,i)*matrix_get(A,k,j);
			matrix_set(R,i,j,R_ij);
			for(int k=0;k<n;k++) {
				a_jk = matrix_get(A,k,j)-matrix_get(A,k,i)*R_ij;
				matrix_set(A,k,j,a_jk);
			}
		}
    }
}
