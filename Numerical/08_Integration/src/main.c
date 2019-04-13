#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<unistd.h>
#include<assert.h>

/*--- TEST FUNCTIONS ---*/


/*--- FUNCTION DECLARATIONS ---*/

double adapt_recur(double f(double x), double f2, double f3, double a, double b, double acc, double eps, int nrec) {
	assert(nrec < 100000);
	double f1 = f(a+1./6*(b-a)), f4 = f(a+5./6*(b-a));
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a), q = (f1+f2+f3+f4)/4*(b-a);
	double tolerance = acc+eps*fabs(Q), error = fabs(Q-q);
	if (error < tolerance) return Q;
	else {
		double Q1 = adapt_recur(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1);
		double Q2 = adapt_recur(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1);
		return Q1+Q2;
	}
}

double adapt_init(double f(double x), double a, double b, double acc, double eps) {
	int nrec = 0;
	double f2 = f(a+1./3*(b-a)), f3 = f(a+2./3*(b-a));
	return adapt_recur(f,f2,f3,a,b,acc,eps,nrec);
}


/*--- MAIN PROGRAM ---*/

int main() {
	
	// A - Recursive adaptive integrator
	
	int calls = 0;
	double a = 0, b = 1, acc = 1e-6, eps = 1e-6, Q;
	double f1(double x) {calls++; return sqrt(x);}
	printf("Integrating sqrt(x) from 0 to 1:\n");
	printf("(acc = eps = 1e-6)\n");
	Q = adapt_init(f1,a,b,acc,eps);
	printf("Q = %g, function calls = %d\n",Q,calls);

	calls = 0;
	double f2(double x) {calls++; return 1/sqrt(x);}
	printf("\nIntegrating 1/sqrt(x) from 0 to 1:\n");
	printf("(acc = eps = 1e-6)\n");
	Q = adapt_init(f2,a,b,acc,eps);
	printf("Q = %g, function calls = %d\n",Q,calls);
	
	calls = 0;
	double f3(double x) {calls++; return log(x)/sqrt(x);}
	printf("\nIntegrating ln(x)/sqrt(x) from 0 to 1:\n");
	printf("(acc = eps = 1e-6)\n");
	Q = adapt_init(f3,a,b,acc,eps);
	printf("Q = %g, function calls = %d\n",Q,calls);
	
	calls = 0, acc = DBL_EPSILON, eps = DBL_EPSILON;
	double f4(double x) {calls++; return 4*sqrt(1-(1-x)*(1-x));}
	printf("\nIntegrating 4*sqrt(1-(1-x)^2) from 0 to 1:\n");
	printf("(acc = eps = DBL_EPSILON)\n");
	Q = adapt_init(f4,a,b,acc,eps);
	printf("Q = %17.16g, function calls = %d\n",Q,calls);
	
	// b - Clenshaw-Curtis variable transformation
	
	
	
	// C - Infinite limits
	
	
	return 0;
}
