#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_linalg.h>
#include"qspline.h"

int binary_search(int n, double *x, double z);

qspline* qspline_alloc(int n, double *x, double *y) {
	qspline* s = malloc(sizeof(qspline));
	s->n = n;
	
	// allocate A x = b
	gsl_matrix *A = gsl_matrix_calloc(3*(n-1),3*(n-1));
	gsl_vector *b = gsl_vector_calloc(3*(n-1));
	gsl_vector *X = gsl_vector_calloc(3*(n-1));
	
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
	
	int j;
	// set 1st derivative equations
	for(int i=2*(n-1);i<3*(n-1)-1;i++) {
		j = 3*(i-2*(n-1));
		gsl_matrix_set( A , i , j+1 , 1                 );
		gsl_matrix_set( A , i , j+2 , 2*x[i-2*(n-1)+1]  );
		gsl_matrix_set( A , i , j+4 , -1                );
		gsl_matrix_set( A , i , j+5 , -2*x[i-2*(n-1)+1] );
	}
	
	// set last equation, derivative at x[1] should be average of derivatives at x[0] and x[2]
	gsl_matrix_set( A , 3*(n-1)-1 , 1 , 1-1.0/2     );
	gsl_matrix_set( A , 3*(n-1)-1 , 2 , 2*x[1]-x[0] );
	gsl_matrix_set( A , 3*(n-1)-1 , 4 , -1.0/2      );
	gsl_matrix_set( A , 3*(n-1)-1 , 5 , -x[2]       );

	// solve linear system
	gsl_linalg_HH_solve(A,b,X);
	
	// allocate x and y data, and a, b, and c coefficients
	s->x = malloc(n*sizeof(double));
	s->y = malloc(n*sizeof(double));
	s->a = malloc((n-1)*sizeof(double));
	s->b = malloc((n-1)*sizeof(double));
	s->c = malloc((n-1)*sizeof(double));
	
	for(int i=0;i<n;i++) {
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

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
	int i = binary_search(s->n,s->x,z);
	return s->a[i] + s->b[i]*z + s->c[i]*z*z;
}

double qspline_deriv(qspline *s, double z) {
	int i = binary_search(s->n,s->x,z);
	return s->b[i] + s->c[i]*2*z;
}

double qspline_integ(qspline *s, double z) {
	int i = binary_search(s->n,s->x,z), j = 0;
	double sum = 0;
	double *x = s->x, *a = s->a, *b = s->b, *c = s->c;
	
	while(j<=i) {
		if(j==i) sum += (a[j]*z    + b[j]*z*z/2       + c[j]*z*z*z/3)
					   -(a[j]*x[j] + b[j]*x[j]*x[j]/2 + c[j]*x[j]*x[j]*x[j]/3);

		else sum += (a[j]*x[j+1] + b[j]*x[j+1]*x[j+1]/2 + c[j]*x[j+1]*x[j+1]*x[j+1]/3)
				   -(a[j]*x[j]   + b[j]*x[j]*x[j]/2     + c[j]*x[j]*x[j]*x[j]/3);
		j++;
	}
	return sum;
}

void qspline_free(qspline *s){
	free(s->x);
	free(s->y);
	free(s->a);
	free(s->b);
	free(s->c);
	free(s);
}
