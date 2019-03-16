#include<stdio.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#define TYPE gsl_multiroot_fdfsolver_hybridsj
#define EPS 1e-6

int tanroot_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double X = *(double*)params;
	double a = gsl_vector_get(x,0);
	gsl_vector_set(f,0,tan(a)-X);
	return GSL_SUCCESS;
}

int tanroot_df (const gsl_vector * x, void * params, gsl_matrix * J)
{
   double a = gsl_vector_get(x,0);
   gsl_matrix_set(J,0,0,1+(tan(a)*tan(a)));
   return GSL_SUCCESS;
}

int tanroot_fdf (const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J)
{
   tanroot_f(x,params,f);
   tanroot_df(x,params,J);
   return GSL_SUCCESS;
}

double myarctan (double x){
	
	// define starting point
	double a0;
	if(x>0) a0 = M_PI/(2+1/x);
	else if(x<0) a0 = -M_PI/(2-1/x);
	else a0 = 0;

	int n = 1;
	gsl_multiroot_function_fdf F;
	F.f      = tanroot_f;
	F.df     = tanroot_df;
	F.fdf    = tanroot_fdf;
	F.n      = n;
	F.params = (void*)&x;

	gsl_multiroot_fdfsolver * S;
	S = gsl_multiroot_fdfsolver_alloc(TYPE,n);
	
	gsl_vector* start = gsl_vector_alloc(n);
	gsl_vector_set(start,0,a0);
	gsl_multiroot_fdfsolver_set(S,&F,start);
	
	int iter = 0, flag;
	
	do{
		iter++;
		gsl_multiroot_fdfsolver_iterate(S);
		flag = gsl_multiroot_test_residual(S->f,EPS);
	}
	while(flag==GSL_CONTINUE);
	
	double a = gsl_vector_get(S->x,0);

	gsl_multiroot_fdfsolver_free(S);
	gsl_vector_free(start);
	
	return a;
}
