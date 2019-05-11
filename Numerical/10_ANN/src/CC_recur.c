#include<stdio.h>
#include<assert.h>
#include<math.h>

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
