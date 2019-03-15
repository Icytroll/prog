#include<stdio.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#define TYPE gsl_multiroot_fsolver_hybrids
#define EPS 1e-12

int tanroot
(const gsl_vector * x, void * params, gsl_vector * f)
{
	double X = *(double*)params;
	double a = gsl_vector_get(x,0);
	gsl_vector_set(f,0,tan(a)-X);
	return GSL_SUCCESS;
}

double myarctan (double x){
	
	int n = 1;
	gsl_multiroot_function F;
	F.f      = tanroot;
	F.n      = n;
	F.params = (void*)&x;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(TYPE,n);

	gsl_vector* start = gsl_vector_alloc(n);
	gsl_vector_set(start,0,0);
	gsl_multiroot_fsolver_set(S,&F,start);
	
	int iter = 0, flag;
	
	do{
		
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);

	}while(flag==GSL_CONTINUE);
	
	double a = gsl_vector_get(S->x,0);

	gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	
	return a;
}
