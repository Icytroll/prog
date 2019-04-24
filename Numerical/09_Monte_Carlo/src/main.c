#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<unistd.h>
#include<assert.h>
#include"vector.h"
#include"matrix.h"

/*--- FUNCTION DECLARATIONS ---*/

double adapt_recur(double f(double x), double f2, double f3, double a, double b, double acc, double eps, int nrec) {
	assert(nrec < 100000);
	assert(~isinf(f2));
	assert(~isinf(f3));
	double f1 = f(a+1./6*(b-a)), f4 = f(a+5./6*(b-a));
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a), q = (f1+f2+f3+f4)/4*(b-a);
	double tol = acc+eps*fabs(Q), error = fabs(Q-q);
	if (error < tol) return Q;
	else {
		double Q1 = adapt_recur(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1);
		double Q2 = adapt_recur(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1);
		return Q1+Q2;
	}
}

double CC_recur(double f(double x), double f2, double f3, double a, double b, double A, double B, double acc, double eps, int nrec) {
	assert(nrec < 100000);
	assert(~isinf(f2));
	assert(~isinf(f3));
	double tet1 = A+1./6*(B-A), tet2 = A+1./3*(B-A), tet3 = A+2./3*(B-A), tet4 = A+5./6*(B-A);
	double f1 = f( a+(b-a)/2. * (cos(tet1) + 1) ), f4 = f( a+(b-a)/2. * (cos(tet4) + 1) );
	double I1 = f1*(b-a)/2.*sin(tet1), I2 = f2*(b-a)/2.*sin(tet2), I3 = f3*(b-a)/2.*sin(tet3), I4 = f4*(b-a)/2.*sin(tet4);
	
	// If one of the integrands are negative or positive infinity, it means we're so close to a divergence point that we don't have enough precision to go deeper. In such a case just return 0, as the subdivision is almost infinitely small so this one value doesn't matter.
	if(isinf(I1)||isinf(I2)||isinf(I3)||isinf(I4)) {
		printf("Encountered a divergence point, returning subdivision Q = 0.\n");
		printf("Continuing ...\n");
		return 0;
	}
	double Q = (2*I1+I2+I3+2*I4)/6*(B-A), q = (I1+I2+I3+I4)/4*(B-A);
	double tol = acc+eps*fabs(Q), error = fabs(Q-q);
	if (error < tol) return Q;
	else {
		double Q1 = CC_recur(f,f1,f2,a,b,A,(A+B)/2,acc/sqrt(2),eps,nrec+1);
		double Q2 = CC_recur(f,f3,f4,a,b,(A+B)/2,B,acc/sqrt(2),eps,nrec+1);
		return Q1+Q2;
	}
}

double adapt_init(double f(double x), double a, double b, double acc, double eps) {
	int nrec = 0;
	double f2 = f(a+1./3*(b-a)), f3 = f(a+2./3*(b-a));
	return adapt_recur(f,f2,f3,a,b,acc,eps,nrec);
}

double CC_init(double f(double x), double a, double b, double acc, double eps) {
	int nrec = 0;
	// transformed limits
	double A = 0, B = M_PI;
	double f2 = f( a+(b-a)/2. * (cos(A+1./3*(B-A)) + 1)), f3 = f( a+(b-a)/2. * (cos(A+2./3*(B-A)) + 1));
	return CC_recur(f,f2,f3,a,b,A,B,acc,eps,nrec);
}



void plainMC(
	double f(vector* x),
	vector* a,
	vector* b,
	int N,
	double* result,
	double* error);


/*--- MAIN PROGRAM ---*/

int main() {
	
	// A - Plain Monte Carlo integration
	
	int n = 3;
	vector* a = vector_alloc(n);
	vector* b = vector_alloc(n);

	double res, err;
	double* result = &res;
	double* error = &err;
	
	vector_set_all(a,0);
	vector_set(b,0,1);
	vector_set(b,1,2);
	vector_set(b,2,3);
	double f1(vector* x) {return 1;};
	int N = 10000;
	plainMC(f1,a,b,N,result,error);
	printf("1 integral:\nResult = %g, error = %g\n",*result,*error);
	
	vector_set_all(a,0);
	vector_set(b,0,1);
	vector_set(b,1,1);
	vector_set(b,2,1);
	double fx(vector* x) {return vector_get(x,0)*vector_get(x,1)*vector_get(x,2);};
	N = 10000000;
	plainMC(fx,a,b,N,result,error);
	printf("x integral:\nResult = %g, error = %g\n",*result,*error);
	
	double fx2(vector* x) {return pow(vector_get(x,0),2)*pow(vector_get(x,1),2)*pow(vector_get(x,2),2);};
	N = 10000000;
	plainMC(fx2,a,b,N,result,error);
	printf("x^2 integral:\nResult = %g, error = %g\n",*result,*error);
	/*
	vector_set_all(a,0);
	vector_set_all(b,M_PI);
	double fcos(vector* x) {return 1./(1-cos(vector_get(x,0))*cos(vector_get(x,1))*cos(vector_get(x,2)));};
	int N = 10000;
	plainMC(fcos,a,b,N,result,error);
	printf("Difficult singular integral:\nResult = %g, error = %g\n",*result,*error);
	*/
	// B - Behaviour of error
	
	// C - Stratified sampling
	
	
	return 0;
}
