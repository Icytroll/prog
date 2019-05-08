#include"matrix.h"
#include"vector.h"
#include<unistd.h>
#include<math.h>

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);
void mv_mult(int transposeA, matrix* A, vector* b, vector* c);
void vector_outer(vector* a, vector* b, matrix* C);

void numerical_gradient(
	double f(vector* x),
	vector* df,
	vector* x,
	double dx) {

	int	n = x->size;
	double fp,fm;
	for(int i=0;i<n;i++) {
		vector_set(x,i,vector_get(x,i)+dx);
		fp = f(x);
		vector_set(x,i,vector_get(x,i)-2*dx);
		fm = f(x);
		vector_set(x,i,vector_get(x,i)+dx);
		vector_set(df,i,(fp-fm)/(2*dx));
	}
}

void broyden_numerical(
	double f(vector* x),
	vector* x,
	double dx,
	double epsilon) {
	
	int n = x->size,maxiter = 10000,iter;
	double alpha = 1e-4, lambda, f_x, f_trial;
	
	// Declare vectors and matrices
	vector* df = vector_alloc(n);
	vector* df_trial = vector_alloc(n);
	vector* y = vector_alloc(n);
	vector* By = vector_alloc(n);
	vector* u = vector_alloc(n);
	matrix* dB = matrix_alloc(n,n);
	matrix* B = matrix_alloc(n,n);
	vector* Dx = vector_alloc(n);
	vector* x_trial = vector_alloc(n);
	
	// Initialize step Dx, inverse Hessian B and gradient df
	vector_set_all(Dx,1);
	matrix_set_identity(B);
	numerical_gradient(f,df,x,dx);
	
	iter = 0;
	while(norm(df) > epsilon && iter<maxiter && norm(Dx)>dx) {
		// Get step size
		vector_mult_all(df,-1);
		mv_mult(0,B,df,Dx);
		vector_mult_all(df,-1);
		// Check if we need to decrease the step size using the Armijo condition
		lambda = 1;
		f_x = f(x);
		vector_add(x,Dx,x_trial);
		f_trial = f(x_trial);
		while(f_trial > f_x+alpha*vector_inner(Dx,df) && lambda > 1.0/pow(2,14)) {
			lambda /= 2;
			vector_mult_all(Dx,0.5);
			vector_add(x,Dx,x_trial);
			f_trial = f(x_trial);
		}
		// If the step is bad, reset the inverse Hessian
		if (lambda == 1.0/pow(2,14)) matrix_set_identity(B);
		// Otherwise update the inverse Hessian
		else {
			numerical_gradient(f,df_trial,x_trial,dx);
			vector_sub(df_trial,df,y);
			mv_mult(0,B,y,By);
			vector_sub(Dx,By,u);
			vector_mult_all(u,1/vector_inner(Dx,y));
			vector_outer(u,Dx,dB);
			matrix_add(B,dB,B);
		}
		// Update the step and get gradient for next iteration
		vector_add(x,Dx,x);
		numerical_gradient(f,df,x,dx);
		iter++;
	}
	printf("Iterations = %d\n",iter);

	// clean up
	vector_free(df);
	vector_free(df_trial);
	vector_free(y);
	vector_free(By);
	vector_free(u);
	matrix_free(dB);
	matrix_free(B);
	vector_free(Dx);
	vector_free(x_trial);
}
