#include<math.h>

double CC_recur(double f(double x), double f2, double f3, double a, double b, double A, double B, double acc, double eps, int nrec);

double CC_init(double f(double x), double a, double b, double acc, double eps) {
	int nrec = 0;
	// transformed limits
	double A = 0, B = M_PI;
	double f2 = f( a+(b-a)/2. * (cos(A+1./3*(B-A)) + 1)), f3 = f( a+(b-a)/2. * (cos(A+2./3*(B-A)) + 1));
	return CC_recur(f,f2,f3,a,b,A,B,acc,eps,nrec);
}
