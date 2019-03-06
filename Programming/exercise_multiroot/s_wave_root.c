#include<stdio.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#define TYPE gsl_multiroot_fsolver_hybrids
#define EPS 1e-12

double fe(double eps, double rmax, FILE* stream);

int master(const gsl_vector* x, void* params, gsl_vector* f) {
	double rmax = *(double*)params;
	double eps = gsl_vector_get(x,0);
	double fval = fe(eps,rmax,NULL);
	gsl_vector_set(f,0,fval);

	return GSL_SUCCESS;
}

double s_wave_root(double eps0, double rmax){
	
	int n = 1;
	gsl_multiroot_function F;
	F.f      = master;
	F.n      = n;
	F.params = (void*)&rmax;

	gsl_multiroot_fsolver * S;
	S = gsl_multiroot_fsolver_alloc(TYPE,n);

	gsl_vector* start = gsl_vector_alloc(n);
	gsl_vector_set(start,0,eps0);
	gsl_multiroot_fsolver_set(S,&F,start);
	
	int iter = 0, flag;
	
	do{
		
		iter++;
		gsl_multiroot_fsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
		
	}while(flag==GSL_CONTINUE);
	
	double f_final = gsl_vector_get(S->x,0);
	fprintf(stderr,"S-wave iterations with rmax = %g : %i\n",rmax,iter);

	gsl_multiroot_fsolver_free(S);
	gsl_vector_free(start);
	
	return f_final;
}
