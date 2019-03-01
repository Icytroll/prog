#include<stdio.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#define TYPE gsl_multiroot_fdfsolver_hybridsj
#define EPS 1e-12

int rosenbrock_f
(const gsl_vector * x, void * params, gsl_vector * f)
{
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);
	gsl_vector_set(f,0,-2*(1-X)-400*X*(Y-X*X));
	gsl_vector_set(f,1,200*(Y-X*X));
	return GSL_SUCCESS;
}

int rosenbrock_df
(const gsl_vector * x, void * params, gsl_matrix * J)
{
	double X = gsl_vector_get(x,0);
	gsl_matrix_set(J,0,0,2+1200*X*X);
	gsl_matrix_set(J,0,1,-400*X);
	gsl_matrix_set(J,1,0,-400*X);
	gsl_matrix_set(J,1,1,200);
	return GSL_SUCCESS;
}

int rosenbrock_fdf
(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J)
{
	double X = gsl_vector_get(x,0);
	double Y = gsl_vector_get(x,1);

	gsl_vector_set(f,0,-2*(1-X)-400*X*(Y-X*X));
	gsl_vector_set(f,1,200*(Y-X*X));
	
	gsl_matrix_set(J,0,0,2+1200*X*X);
	gsl_matrix_set(J,0,1,-400*X);
	gsl_matrix_set(J,1,0,-400*X);
	gsl_matrix_set(J,1,1,200);
	
	return GSL_SUCCESS;
}


double rosenbrock(double x, double y) {
	return (1-x)*(1-x) + 100*(y-x*x)*(y-x*x);
}

int rosen_extremum (double x0, double y0){
	
	int n = 2;
	gsl_multiroot_function_fdf F;
	F.f      = rosenbrock_f;
	F.df     = rosenbrock_df;
	F.fdf    = rosenbrock_fdf;
	F.n      = n;
	F.params = (void*)NULL;

	gsl_multiroot_fdfsolver * S;
	S = gsl_multiroot_fdfsolver_alloc(TYPE,n);

	gsl_vector* start = gsl_vector_alloc(n);
	gsl_vector_set(start,0,x0);
	gsl_vector_set(start,1,y0);
	gsl_multiroot_fdfsolver_set(S,&F,start);
	
	int iter = 0, flag;
	double x, y;
	
	x = gsl_vector_get(S->x,0);
	y = gsl_vector_get(S->x,1);

	fprintf(stdout,"%g %g %g\n",x,y,rosenbrock(x,y));
	do{
		
		iter++;
		gsl_multiroot_fdfsolver_iterate(S);
		flag=gsl_multiroot_test_residual(S->f,EPS);
		
		x = gsl_vector_get(S->x,0);
		y = gsl_vector_get(S->x,1);
		fprintf(stdout,"%g %g %g\n",x,y,rosenbrock(x,y));

	}while(flag==GSL_CONTINUE);
	
	fprintf(stderr,"Iterations = %i\n",iter);

	gsl_multiroot_fdfsolver_free(S);
	gsl_vector_free(start);
	
	return flag;
}
