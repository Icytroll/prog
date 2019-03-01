#include<#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#define TYPE gsl_multiroot_fsolver_hybrids
#define EPS 1e-12

int root_equation
(const gsl_vector * x, void * params, gsl_vector * f)
{
	double z = *(double*)params;
	double r = gsl_vector_get(x,0);
	double mismatch = r*r-z;
	gsl_vector_set(f,0,mismatch);
return GSL_SUCCESS;
}

double rosenblock (double x, double y){
	gsl_multiroot_function F;
	F.f=root_equation;
	F.n=1;
	F.params=(void*)&z;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(TYPE,F.n);

	gsl_vector* start = gsl_vector_alloc(F.n);
	gsl_vector_set(start,0,z/2);
	gsl_multiroot_fsolver_set(S,&F,start);

	int flag;
	do{
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
	}while(flag==GSL_CONTINUE);

	double result=gsl_vector_get(S->x,0);
	gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	return result;
}
