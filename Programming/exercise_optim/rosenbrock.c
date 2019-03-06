#include<gsl/gsl_multimin.h>
#include<gsl/gsl_vector.h>

double fval(double x, double y) {
	double res = (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
	return res;
}

double master(const gsl_vector *X, void *params){
	double x = gsl_vector_get(X,0);
	double y = gsl_vector_get(X,1);
	double res = fval(x,y);
	return res;
}


void rosenbrock(double x0, double y0) {

	int dim = 2;
	gsl_multimin_function F;
	F.f = master;
	F.n = dim;
	F.params=NULL;

	gsl_multimin_fminimizer *M;
	#define TYPE gsl_multimin_fminimizer_nmsimplex2
	M = gsl_multimin_fminimizer_alloc(TYPE,dim);
	gsl_vector* start = gsl_vector_alloc(dim);
	gsl_vector* step  = gsl_vector_alloc(dim);
	gsl_vector_set(start,0,x0);
	gsl_vector_set(start,1,y0);
	gsl_vector_set(step,0,0.1);
	gsl_vector_set(step,1,0.1);

	gsl_multimin_fminimizer_set(M,&F,start,step);

	int iter=0,status;
	double size;
	printf("#Iter:  x:      y:      fval:   size:\n");
	printf ("%7.5d %7.5g %7.5g %7.5g %7.5s\n",iter,x0,y0,fval(x0,y0),"NaN");
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(M);
		if (status) break;

		size = gsl_multimin_fminimizer_size (M);
		status = gsl_multimin_test_size (size, 1e-4);

		printf ("%7.5d %7.5g %7.5g %7.5g %7.5g\n",
              iter,
              gsl_vector_get (M->x, 0),
              gsl_vector_get (M->x, 1),
              M->fval, size);
    }
	while (status == GSL_CONTINUE && iter < 1000);

	gsl_vector_free(start);
	gsl_vector_free(step);
	gsl_multimin_fminimizer_free(M);
}
