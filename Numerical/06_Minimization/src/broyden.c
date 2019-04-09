#include"matrix.h"
#include"vector.h"
#include<math.h>

void qr_gs_decomp(matrix* Q, matrix* R);
void qr_gs_solve(matrix* Q, matrix* R, vector* b, vector* x);
void mv_mult(int transposeA, matrix* A, vector* b, vector* c);
void vector_outer(vector* a, vector* b, matrix* C);


void broyden(
	double f(vector* x),
	void set_df(vector* x, vector* df),
	vector* x,
	double epsilon) {
	
	int n = x->size,maxiter = 1000,iter;
	double alpha = 1e-4, lambda, f_x, f_trial;
	vector* df = vector_alloc(n);
	vector* df_trial = vector_alloc(n);
	vector* y = vector_alloc(n);
	vector* By = vector_alloc(n);
	vector* u = vector_alloc(n);
	matrix* dB = matrix_alloc(n,n);
	matrix* B = matrix_alloc(n,n);
	vector* Dx = vector_alloc(n);
	vector* x_trial = vector_alloc(n);
	
	vector_set_all(Dx,1);
	matrix_set_identity(B);
	
	set_df(x,df);
	iter = 0;
	while(norm(df) > epsilon && iter<maxiter) {
		vector_mult_all(df,-1);
		mv_mult(0,B,df,Dx);
		vector_mult_all(df,-1);
		
		lambda = 1;
		f_x = f(x);
		vector_add(x,Dx,x_trial);
		f_trial = f(x_trial);
		//vector_print(x_trial,"x_trial =",stdout);
		//printf("f_trial = %g, f_x+alpha*s'*df = %g\n",f_trial,f_x+alpha*vector_inner(Dx,df));
		while(f_trial > f_x+alpha*vector_inner(Dx,df) && lambda > 1.0/pow(2,14)) {
			lambda /= 2;
			vector_mult_all(Dx,0.5);
			vector_add(x,Dx,x_trial);
			f_trial = f(x_trial);
			//vector_print(x_trial,"x_trial =",stdout);
			//printf("f_trial = %g, f_x+alpha*s'*df = %g\n",f_trial,f_x+alpha*vector_inner(Dx,df));
			//printf("lambda = %g\n",lambda);
		}
		if (lambda == 1.0/pow(2,14)) matrix_set_identity(B);
		else {
			set_df(x_trial,df_trial);
			//vector_print(df_trial,"df_trial =",stdout);
			vector_sub(df_trial,df,y);
			//vector_print(y,"y =",stdout);
			mv_mult(0,B,y,By);
			//vector_print(By,"B*y =",stdout);
			vector_sub(Dx,By,u);
			//vector_print(u,"u =",stdout);
			vector_mult_all(u,1/vector_inner(Dx,y));
			//vector_print(u,"c =",stdout);
			vector_outer(u,Dx,dB);
			//matrix_print(dB,"dB =",stdout);
			matrix_add(B,dB,B);
			//matrix_print(B,"B =",stdout);
		}
		
		vector_add(x,Dx,x);
		//vector_print(x,"x =",stdout);
		set_df(x,df);
		//vector_print(df,"df =",stdout);
		iter++;
	}
	fprintf(stderr,"%10d\n",iter);

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
