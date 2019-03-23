#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"matrix.h"
#include"vector.h"

void qr_gs_inverse(matrix* Q, matrix* R, matrix* B) {
	
	// Ax_i = e_i --> QRx_i = e_i --> Rx_i = QTe_i --> Rx_i = QT_i
	// With e_i being the ith coloumn of an identity matrix, x_i
	// will be the ith coloumn of the inverse of A.

	int n = Q->size1;
	matrix* QT = matrix_transpose(Q);
	
	double s;
	for(int i=0;i<n;i++) {
		for(int j=n-1;j>=0;j--) {
			s = matrix_get(QT,j,i);
			for(int k=j+1;k<n;k++) {
				s -= matrix_get(R,j,k)*matrix_get(B,k,i);
			}
			matrix_set(B,j,i,s/matrix_get(R,j,j));
		}
	}
	matrix_free(QT);
}
