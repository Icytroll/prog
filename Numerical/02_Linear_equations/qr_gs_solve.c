#include"matrix.h"
#include"vector.h"

void qr_gs_solve(matrix* QT, matrix* R, vector* b, vector* x) {
	int n = QT->size1;
	
	// Apply QT to b
	double sum;
    for(int i=0;i<n;i++) {
        sum = 0;
        for(int j=0;j<n;j++) {
			sum += matrix_get(QT,i,j)*vector_get(b,j);
		}
		vector_set(x,i,sum);
    }
	// Perform back-substitution
	for(int i=n-1;i>=0;i--) {
		double s = vector_get(x,i);
		for(int j=i+1;j<n;j++)
			s -= matrix_get(R,i,j)*vector_get(x,j);
		vector_set(x,i,s/matrix_get(R,i,i));
	}
}
