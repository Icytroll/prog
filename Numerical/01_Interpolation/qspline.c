#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include"qspline.h"

int binary_search(int n, double *x, double z);

qspline* qspline_alloc(int n, double *x, double *y) {
	printf("Alloc 1\n");
	qspline* s = malloc(sizeof(qspline));
	s->n = n;
	
	printf("Alloc 2\n");
	// allocate A x = b
	gsl_matrix *A = gsl_matrix_calloc(3*(n-1),3*(n-1));
	gsl_vector *b = gsl_vector_calloc(3*(n-1));
	gsl_vector *X = gsl_vector_calloc(3*(n-1));
	
	printf("Alloc 3\n");
	// set interpolation (node) equations
	for(int i=0;i<n-1;i++) {
		gsl_matrix_set( A , 2*i   , 3*i   , 1             );
		gsl_matrix_set( A , 2*i   , 3*i+1 , x[i]          );
		gsl_matrix_set( A , 2*i   , 3*i+2 , x[i]*x[i]     );
		gsl_matrix_set( A , 2*i+1 , 3*i   , 1             );
		gsl_matrix_set( A , 2*i+1 , 3*i+1 , x[i+1]        );
		gsl_matrix_set( A , 2*i+1 , 3*i+2 , x[i+1]*x[i+1] );
		
		gsl_vector_set( b , 2*i   , y[i]   );
		gsl_vector_set( b , 2*i+1 , y[i+1] );
	}
	
	printf("Alloc 4\n");
	int j;
	// set 1st derivative equations
	for(int i=2*(n-1);i<3*(n-1)-1;i++) {
		j = 3*(i-2*(n-1));
		printf("i=%i j=%i\n",i,j);
		gsl_matrix_set( A , i , j+1 , 1       );
		gsl_matrix_set( A , i , j+2 , x[i-2*(n-1)+1]  );
		gsl_matrix_set( A , i , j+4 , -1      );
		gsl_matrix_set( A , i , j+5 , -x[i-2*(n-1)+1] );
	}
	
	printf("Alloc 5\n");
	// set last equation, c1 = 0 (first spline is assumed linear)
	gsl_matrix_set( A , 3*(n-1)-1 , 2 , 1);


	for(int i=0;i<3*(n-1);i++) {
		for(int j=0;j<3*(n-1);j++) {
			printf("%g ",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}

	printf("Alloc 6\n");
	// solve linear system
	gsl_linalg_HH_solve(A,b,X);
	
	/*
	int r;
	gsl_permutation * p = gsl_permutation_alloc (n-1);
	gsl_linalg_LU_decomp (A, p, &r);
	gsl_linalg_LU_solve  (A, p, b, X);
	*/
	
	printf("Alloc 7\n");
	// allocate a, b, and c coefficients
	s->a = malloc((n-1)*sizeof(double));
	s->b = malloc((n-1)*sizeof(double));
	s->c = malloc((n-1)*sizeof(double));
	
	// copy coefficients to our structure
	for(int i=0;i<n-1;i++) {
		s->a[i] = gsl_vector_get(X,3*i);
		s->b[i] = gsl_vector_get(X,3*i+1);
		s->c[i] = gsl_vector_get(X,3*i+2);
	}

	// free up memory
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(X);
//	gsl_permutation_free(p);
		
	return s;
}


double qspline_eval(qspline *s, double z) {
	printf("Eval 1\n");
	int n = s->n;
	double *x = s->x;
	int i = binary_search(n,s->x,z);
	printf("Eval 2\n");
	return s->a[i] + s->b[i]*z + s->c[i]*z*z;
}

double qspline_deriv(qspline *s, double z) {
	int i = binary_search(s->n,s->x,z);
	return s->b[i] + s->c[i]*z;
}

double qspline_integ(qspline *s, double z) {
	int i = binary_search(s->n,s->x,z), j = 0;
	double sum = 0;
	double *x = s->x, *a = s->a, *b = s->b, *c = s->c;
	
	while(j<=i) {
		if(j==i) sum += (a[j]*z    + b[j]*z*z/2       + c[j]*z*z*z/6)
					   -(a[j]*x[j] + b[j]*x[j]*x[j]/2 + c[j]*x[j]*x[j]*x[j]/6);

		else sum += (a[j]*x[j+1] + b[j]*x[j+1]*x[j+1]/2 + c[j]*x[j+1]*x[j+1]*x[j+1]/6)
				   -(a[j]*x[j]   + b[j]*x[j]*x[j]/2     + c[j]*x[j]*x[j]*x[j]/6);
		j++;
	}
	return sum;
}

void qspline_free(qspline *s){
	free(s->a);
	free(s->b);
	free(s->c);
	free(s);
}
