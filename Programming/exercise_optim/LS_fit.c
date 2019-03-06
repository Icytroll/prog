#include<stdio.h>
#include<gsl/gsl_multimin.h>
#include<gsl/gsl_vector.h>
#include<math.h>

double Fval(int n, double A, double T, double B, double* t, double* y, double* e) {
	double sum = 0;
	#define f(t) A*exp(-(t)/T) + B
	for(int i=0;i<n;i++) sum += pow( (f(t[i]) - y[i])/e[i] ,2);
	return sum;
}

struct expdata {int n; double *t,*y,*e;};

double fit_function (const gsl_vector *x, void *params) {
	double  A = gsl_vector_get(x,0);
	double  T = gsl_vector_get(x,1);
	double  B = gsl_vector_get(x,2);
	struct expdata *p = (struct expdata*) params;
	int     n = p->n;
	double *t = p->t;
	double *y = p->y;
	double *e = p->e;
	return Fval(n,A,T,B,t,y,e);
}


void LS_fit(int n, double* t, double* y, double* e) {
	
	int dim=3;
	struct expdata data;
	data.n=n;
	data.t=t;
	data.y=y;
	data.e=e;
	gsl_multimin_function F;
	F.f = fit_function;
	F.n = dim;
	F.params=(void*)&data;

	gsl_multimin_fminimizer *M;
	#define TYPE gsl_multimin_fminimizer_nmsimplex2
	M = gsl_multimin_fminimizer_alloc(TYPE,dim);
	gsl_vector* start=gsl_vector_alloc(dim);
	gsl_vector* step=gsl_vector_alloc(dim);
	gsl_vector_set(start,0,1);
	gsl_vector_set(start,1,1);
	gsl_vector_set(start,2,1);
	gsl_vector_set(step,0,0.5);
	gsl_vector_set(step,1,0.5);
	gsl_vector_set(step,2,0.5);

	gsl_multimin_fminimizer_set(M,&F,start,step);

	int iter=0,status;
	double size;
	fprintf(stderr,"#Iter:  A:      T:      B:      fval:   size:\n");
	fprintf(stderr,"%7.5d %7.5d %7.5d %7.5d %7.5g %7.5s\n",iter,1,1,1,Fval(n,1,1,1,t,y,e),"NaN");
	
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(M);
		if (status) break;

		size = gsl_multimin_fminimizer_size (M);
		status = gsl_multimin_test_size (size, 1e-4);

		fprintf(stderr,"%7.5d %7.5g %7.5g %7.5g %7.5g %7.5g\n",
              iter,
              gsl_vector_get (M->x, 0),
              gsl_vector_get (M->x, 1),
			  gsl_vector_get (M->x, 2),
              M->fval, size);
    }
	while (status == GSL_CONTINUE && iter < 1000);
	
	double A = gsl_vector_get(M->x,0);
	double T = gsl_vector_get(M->x,1);
	double B = gsl_vector_get(M->x,2);
	fprintf(stderr,"\n\n#t:       f:\n");
	for(double i=0;i<=10;i+=0.1)
		fprintf(stderr,"%7.5g %7.5g\n",i,A*exp(-i/T)+B);

	gsl_vector_free(start);
	gsl_vector_free(step);
	gsl_multimin_fminimizer_free(M);	
}
