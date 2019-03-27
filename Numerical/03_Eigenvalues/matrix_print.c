#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>


void matrix_print(gsl_matrix* A, char* s, FILE* stream) {
    fprintf(stream,"%s\n",s);
    int n = A->size1;
    int m = A->size2;
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            fprintf(stream,"%7.3f ",gsl_matrix_get(A,i,j));
        }
        fprintf(stream,"\n");
    }
    fprintf(stream,"\n");
}

