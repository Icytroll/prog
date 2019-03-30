#include"matrix.h"
#include"vector.h"

void mv_mult(int transposeA, matrix* A, vector* b, vector* c) {
    
	int n, m;
	double sum;
    if (transposeA) {
		n = A->size2;
		m = A->size1;
    	matrix* AT = matrix_alloc(m,n);
		matrix_transpose(A,AT);
		for(int i=0;i<n;i++) {
        	sum = 0;
        	for(int j=0;j<m;j++) {
            	sum += matrix_get(AT,i,j)*vector_get(b,j);
			}
        	vector_set(c,i,sum);
    	}
		matrix_free(AT);
	}

	else {
		n = A->size1;
		m = A->size2;
    	for(int i=0;i<n;i++) {
        sum = 0;
        	for(int j=0;j<m;j++) {
            	sum += matrix_get(A,i,j)*vector_get(b,j);
			}
        	vector_set(c,i,sum);
    	}
	}
}
