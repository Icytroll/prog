#include"matrix.h"
#include"vector.h"
#include<math.h>

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);

void numeric_jacobian(
	void f(vector* x, vector* fx),
	matrix* J,
	vector* fx,
	vector* x,
	double dx,
	int* fcptr) {

	int	n = x->size;
	double fp,fm;
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			vector_set(x,j,vector_get(x,j)+dx);
			f(x,fx); (*fcptr)++;
			fp = vector_get(fx,i);
			vector_set(x,j,vector_get(x,j)-2*dx);
			f(x,fx); (*fcptr)++;
			fm = vector_get(fx,i);
			vector_set(x,j,vector_get(x,j)+dx);
			matrix_set(J,i,j,(fp-fm)/(2*dx));
		}
	}
}

void newton(void f(vector* x, vector* fx),
			vector* x,
			double dx,
			double epsilon) {
	
	int n = x->size,maxiter = 1000,iter,fcount = 0;
	int* fcptr = &fcount;
	double lambda, norm_x, norm_trial;
	vector* fx = vector_alloc(n);
	matrix* J = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	vector* Dx = vector_alloc(n);
	vector* x_trial = vector_alloc(n);

	vector_set_all(Dx,1);
	
	f(x,fx); (*fcptr)++;
	iter = 0;
	while(norm(fx)>epsilon && iter<maxiter && norm(Dx)>dx) {
		numeric_jacobian(f,J,fx,x,dx,fcptr);
		qr_gs_decomp(J,R);
		vector_mult_all(fx,-1);
		qr_gs_solve(J,R,fx,Dx);
		vector_mult_all(fx,-1);
		
		lambda = 1;
		vector_add(x,Dx,x_trial);
		norm_x = norm(fx);
		f(x_trial,fx); (*fcptr)++;
		norm_trial = norm(fx);
		while(norm_trial > (1-lambda/2)*norm_x && lambda > 1.0/64) {
			lambda /= 2;
			vector_mult_all(Dx,0.5);
			vector_add(x,Dx,x_trial);
			f(x_trial,fx); (*fcptr)++;
			norm_trial = norm(fx);
		}
		vector_add(x,Dx,x);
		f(x,fx); (*fcptr)++;
		iter++;
	}
	fprintf(stderr,"%10d %10d\n",iter,*fcptr);

	vector_free(fx);
	matrix_free(J);
	matrix_free(R);
	vector_free(Dx);
	vector_free(x_trial);
}
