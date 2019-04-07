#include<stdio.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#include<unistd.h>
#define TYPE gsl_multiroot_fdfsolver_hybridsj

void gsl_rootfinder(
	int f(const gsl_vector* x, void* params, gsl_vector* f),
	int df(const gsl_vector* x, void* params, gsl_matrix* J),
	int fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J),
	gsl_vector* x,
	double epsilon) {
	
	int n = x->size;
	gsl_multiroot_function_fdf F;
	F.f      = f;
	F.df     = df;
	F.fdf    = fdf;
	F.n      = n;
	F.params = (void*)NULL;
	
	gsl_multiroot_fdfsolver* S;
	S = gsl_multiroot_fdfsolver_alloc(TYPE,n);
	gsl_multiroot_fdfsolver_set(S,&F,x);
	
	int iter = 0, flag;
	do{
		iter++;
		gsl_multiroot_fdfsolver_iterate(S);
		flag = gsl_multiroot_test_residual(S->f,epsilon);
	}
	while(flag==GSL_CONTINUE);
	
	fprintf(stderr,"%10d\n",iter);
	gsl_multiroot_fdfsolver_free(S);
	
}
