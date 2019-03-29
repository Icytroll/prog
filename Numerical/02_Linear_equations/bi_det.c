#include"matrix.h"

double bi_det(matrix* B) {
	double det = 1;
	int n = B->size1;
	for(int i=0;i<n;i++) det *= matrix_get(B,i,i);
	return det;
}
