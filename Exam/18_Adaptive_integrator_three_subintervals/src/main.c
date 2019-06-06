#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<unistd.h>
#include<assert.h>

/*

Jens S. K. Jensen - student number: 201209587
Exam question -> mod(87,23) = 18:
	Adaptive integrator with subdivision into three subintervals.

1D adaptive integrator using 3 uniform subdivisions and the following points:
x_i = 1/6, 3/6, 5/6
w_i = 3/8, 2/8, 3/8\n");
For each subdivision the points x_i of the current subdivision will be the new midpoints, meaning a subdivision at the first point x = 1/6 will have the new points:
x_i = 1/18, 3/18, 5/18
meaning we only reuse 1 point per subdivision. We're using an open integration scheme to avoid singularities at the integration limits.

*/


/* -----=====¤¤¤¤¤ FUNCTION DECLARATIONS ¤¤¤¤¤=====----- */

double adapt_recur(double f(double x), double f2, double a, double b, double acc, double eps, int nrec, int plot) {
	assert(nrec < 100000); // stop if we're going too deep in a single subdivision
	assert(~isinf(f2));    // stop if a sample point is on top of a singularity
	double dx = b-a;
	double f1 = f(a+1./6*dx), f3 = f(a+5./6*dx);
	double Q = (3*f1+2*f2+3*f3)/8.*dx, q = (f1+f2+f3)/3.*dx;
	double tol = acc+eps*fabs(Q), error = fabs(Q-q);
	if (plot) fprintf(stderr,"%g %g %g %g\n%g %g %g %g\n\n",a,0.,b,0.,a,f(a),b,f(b)); // used to plot subdivisions
	if (error < tol) return Q;
	else {
		double Q1 = adapt_recur(f,f1, a         ,a+1.*dx/3 ,acc/sqrt(2),eps,nrec+1,plot);
		double Q2 = adapt_recur(f,f2, a+1.*dx/3 ,a+2.*dx/3 ,acc/sqrt(2),eps,nrec+1,plot);
		double Q3 = adapt_recur(f,f3, a+2.*dx/3 ,b         ,acc/sqrt(2),eps,nrec+1,plot);
		return Q1+Q2+Q3;
	}
}

double adapt_init(double f(double x), double a, double b, double acc, double eps, int plot) {
	int nrec = 0;
	double f2 = f(a+3./6*(b-a));
	return adapt_recur(f,f2,a,b,acc,eps,nrec,plot);
}


/* -----=====¤¤¤¤¤ MAIN PROGRAM ¤¤¤¤¤=====----- */

int main() {
	
	// Test a few trial functions
	
	printf("I = exact solution\n");
	printf("Q = numerical approximation\n");
	
	int calls = 0, plot = 0;
	double Q, a = 0, b = 1, acc = 1e-6, eps = 1e-6;
	
	double f1(double x) {calls++; return x;}
	printf("\nIntegrating f(x) = x from 0 to 1 (acc = eps = 1e-6) ...\n");
	Q = adapt_init(f1,a,b,acc,eps,plot);
	printf("I = 0.5\n");
	printf("Q = %.16g, function calls = %d\n",Q,calls);
	
	calls = 0;
	double f2(double x) {calls++; return sqrt(x);}
	printf("\nIntegrating f(x) = sqrt(x) from 0 to 1 (acc = eps = 1e-6) ...\n");
	Q = adapt_init(f2,a,b,acc,eps,plot);
	printf("I = 0.666 ...\n");
	printf("Q = %.16g, function calls = %d\n",Q,calls);

	calls = 0;
	double f3(double x) {calls++; return 1/sqrt(x);}
	printf("\nIntegrating f(x) = 1/sqrt(x) from 0 to 1 (acc = eps = 1e-6) ...\n");
	Q = adapt_init(f3,a,b,acc,eps,plot);
	printf("I = 2\n");
	printf("Q = %.16g, function calls = %d\n",Q,calls);
	
	calls = 0;
	double f4(double x) {calls++; return log(x)/sqrt(x);}
	printf("\nIntegrating f(x) = ln(x)/sqrt(x) from 0 to 1 (acc = eps = 1e-6) ...\n");
	Q = adapt_init(f4,a,b,acc,eps,plot);
	printf("I = -4\n");
	printf("Q = %.16g, function calls = %d\n",Q,calls);

	calls = 0, acc = DBL_EPSILON, eps = DBL_EPSILON;
	double f5(double x) {calls++; return 4*sqrt(1-(1-x)*(1-x));}
	printf("\nIntegrating f(x) = 4*sqrt(1-(1-x)^2) from 0 to 1 (acc = eps = DBL_EPSILON) ...\n");
	Q = adapt_init(f5,a,b,acc,eps,plot);
	printf("I = %.16g ...\n",M_PI);
	printf("Q = %.16g, function calls = %d\n",Q,calls);
	
	calls = 0, a = 0, b = 2*M_PI, acc = 1e-6, eps = 1e-6, plot = 1;
	double f6(double x) {calls++; return sin(x)/x;}
	printf("\nIntegrating f(x) = sin(x)/x from 0 to 2*pi (acc = eps = 1e-6) ...\n");
	Q = adapt_init(f6,a,b,acc,eps,plot);
	printf("I = 1.418151576132628 ...\n");
	printf("Q = %.16g, function calls = %d\n",Q,calls);
	
	return 0;
}
