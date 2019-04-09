#include"matrix.h"
#include"vector.h"
#include<math.h>

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);

void newton(
	double f(vector* x),
	void set_dfH(vector* x, vector* df, matrix* H),
	vector* x,
	double epsilon) {
	
	int n = x->size,maxiter = 1000,iter;
	double alpha = 1e-4, lambda, f_x, f_trial;
	vector* df = vector_alloc(n);
	matrix* H = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	vector* Dx = vector_alloc(n);
	vector* x_trial = vector_alloc(n);

	vector_set_all(Dx,1);
	
	set_dfH(x,df,H);
	iter = 0;
	while(norm(df) > epsilon && iter<maxiter) {
		qr_gs_decomp(H,R);
		vector_mult_all(df,-1);
		qr_gs_solve(H,R,df,Dx);
		vector_mult_all(df,-1);
		
		lambda = 1;
		f_x = f(x);
		vector_add(x,Dx,x_trial);
		f_trial = f(x_trial);
		while(f_trial > f_x+alpha*vector_inner(Dx,df) && lambda > 1.0/64) {
			lambda /= 2;
			vector_mult_all(Dx,0.5);
			vector_add(x,Dx,x_trial);
			f_trial = f(x_trial);
		}
		vector_add(x,Dx,x);
		set_dfH(x,df,H);
		iter++;
	}
	fprintf(stderr,"%10d\n",iter);

	vector_free(df);
	matrix_free(H);
	matrix_free(R);
	vector_free(Dx);
	vector_free(x_trial);
}
