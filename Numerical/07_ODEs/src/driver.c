#include"vector.h"
#include"matrix.h"
#include<math.h>

int driver(
    double* t,
    double b,
    double* h,
    vector* yt,
    double acc, double eps,
    void stepper(
        double t, double h, vector* yt,
        void f(double t, vector* y, vector* dydt, void* params),
        void* params, vector* yth, vector* err),
    void f(double t, vector* y, vector* dydt, void* params),
	void* params,
	matrix* data) {
	
	int n = yt->size, iter = 0, maxiter = data->size1;
	double maxerr;
	vector* yth = vector_alloc(n);
	vector* err = vector_alloc(n);

	for(int i=0;i<n;i++) matrix_set(data,iter,i+1,vector_get(yt,i));
	matrix_set(data,iter,0,*t);

	while (*t<b && iter<maxiter) {
		stepper(*t,*h,yt,f,params,yth,err);
		maxerr = vector_absmax(err);

		if (maxerr < acc || maxerr < eps*maxerr) {
			iter++;
			*t = *t+*h;
			for(int i=0;i<n;i++) {
				vector_set(yt,i,vector_get(yth,i));
				matrix_set(data,iter,i+1,vector_get(yt,i));
			}
			matrix_set(data,iter,0,*t);
			if (maxerr == 0) (*h) = 2*(*h);
			else (*h) = (*h)*pow(acc/maxerr,0.25);
			
		}
		else {
			(*h) = (*h)*pow(acc/maxerr,0.25);
		}
	}
	printf("Iterations: %i\n",iter);
	vector_free(yth);
	vector_free(err);
	return iter;
}
