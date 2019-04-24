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
	
	int n = 3, N;
	vector* a = vector_alloc(n);
	vector* b = vector_alloc(n);
	double res, err;
	double* result = &res;
	double* error = &err;
	
	
	printf("Performing plain Monte Carlo integration:\nN = 10000000\n");
	vector_set_all(a,0);
	vector_set_all(b,1);
	double fx2(vector* x) {return pow(vector_get(x,0),2)*pow(vector_get(x,1),2)*pow(vector_get(x,2),2);};
	N = 1000000;
	plainMC(fx2,a,b,N,result,error);
	printf("\nx^2*y^2*z^2 from x = 0:1, y = 0:1, z = 0:1\nShould be: 1/27 = 0.037037\nResult = %g, error = %g\n",*result,*error);	

	vector_set(b,0,M_PI);
	vector_set(b,1,1);
	vector_set(b,2,M_PI);
	double fxyz(vector* x) {
		double x1 = vector_get(x,0);
		double x2 = vector_get(x,1);
		double x3 = vector_get(x,2);
		return sin(x1)*sin(x1)+x2*sin(x3);
	};
	N = 1000000;
	plainMC(fxyz,a,b,N,result,error);
	printf("\nsin(x)^2+y*sin(z) from x = 0:pi, y = 0:1, z = 0:pi\nShould be: 0.5*pi*(2+pi) = 8.07639\nResult = %g, error = %g\n",*result,*error);
	
	vector_set_all(b,M_PI);
	double fcos(vector* x) {
		double x1 = vector_get(x,0);
		double x2 = vector_get(x,1);
		double x3 = vector_get(x,2);
		return 1/(1-cos(x1)*cos(x2)*cos(x3))/M_PI/M_PI/M_PI;
	};
	N = 1000000;
	plainMC(fcos,a,b,N,result,error);
	printf("\n1/((1-cos(x)*cos(y)*cos(z))*pi^3) from x = 0:pi, y = 0:pi, z = 0:pi\nShould be: G(1/4)^4/(4*pi^3) = 1.39320392\nResult = %g, error = %g\n",*result,*error);
	
	// B - Behaviour of error
	
	int Nlower = 10, Nupper = 100000000;
	for(int i=Nlower;i<=Nupper;i*=Nlower) {
		vector_set_all(a,0);
		vector_set_all(b,1);
		plainMC(fx2,a,b,i,result,error);
		fprintf(stderr,"%d %g",i,*error);
		
		vector_set(b,0,M_PI);
		vector_set(b,2,M_PI);
		plainMC(fxyz,a,b,i,result,error);
		fprintf(stderr," %g",*error);
		
		vector_set(b,1,M_PI);
		plainMC(fcos,a,b,i,result,error);
		fprintf(stderr," %g\n",*error);
	}
		
	// C - Stratified sampling
	
	
	return 0;
}
