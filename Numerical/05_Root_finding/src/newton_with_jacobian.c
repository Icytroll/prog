#include"matrix.h"
#include"vector.h"
#include<math.h>

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);

void newton_with_jacobian(void f(vector* x, vector* fx),void df(vector* x, matrix* J),vector* x,double epsilon) {
	
	int n = x->size,maxiter = 1000,iter;
	double lambda, norm_x, norm_trial;
	vector* fx = vector_alloc(n);
	matrix* J = matrix_alloc(n,n);
	matrix* R = matrix_alloc(n,n);
	vector* Dx = vector_alloc(n);
	vector* x_trial = vector_alloc(n);
	
	f(x,fx);
	iter = 0;
	while(norm(fx)>epsilon && iter<maxiter) {
		df(x,J);
		qr_gs_decomp(J,R);
		vector_mult_all(fx,-1);
		qr_gs_solve(J,R,fx,Dx);
		vector_mult_all(fx,-1);
		
		lambda = 1;
		vector_add(x,Dx,x_trial);
		norm_x = norm(fx);
		f(x_trial,fx);
		norm_trial = norm(fx);
		while(norm_trial > (1-lambda/2)*norm_x && lambda > 1.0/64) {
			lambda /= 2;
			vector_mult_all(Dx,0.5);
			vector_add(x,Dx,x_trial);
			f(x_trial,fx);
			norm_trial = norm(fx);
		}
		vector_add(x,Dx,x);
		f(x,fx);
		iter++;
	}
	printf("Solution found after %d iterations.\n",iter);
	
	vector_free(fx);
	matrix_free(J);
	matrix_free(R);
	vector_free(Dx);
	vector_free(x_trial);
}
