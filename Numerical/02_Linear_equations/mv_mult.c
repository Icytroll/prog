#include"matrix.h"
#include"vector.h"

void mv_mult(matrix* A, vector* b, vector* c) {
    int n = A->size1, m = A->size2;
    double sum;
    for(int i=0;i<n;i++) {
        sum = 0;
        for(int j=0;j<m;j++)
            sum += matrix_get(A,i,j)*vector_get(b,j);
        vector_set(c,i,sum);
    }
}
