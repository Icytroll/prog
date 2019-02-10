#include<stdio.h>
#include<math.h>
#include<complex.h>

int main() {

	// Gamma function
	int n = 5;
	printf("Gamma(%i) = %g\n",n,tgamma(n));

	// Bessel function
	double m = 0.5;
	printf("J1(%g) = %g\n",m,j1(m));

	// Squareroot of -2
	complex l = csqrt(-2);
	printf("Squareroot of -2 = %g + %gi\n",creal(l),cimag(l));
	
	// e^i
	complex o = cexp(I);
	printf("e^i = %g + %gi\n",creal(o),cimag(o));
	
	// e^(i*pi)
	complex p = cexp(I*M_PI);
	printf("e^(i*pi) = %g + %gi\n",creal(p),cimag(p));

	// i^e
	complex q = cpowf(I,M_E);
	printf("i^e = %g + %gi\n",creal(q),cimag(q));
	

	
	// Significant digits for float, double and long double
	float r = 0.1111111111111111111111111;
	double s = 0.1111111111111111111111111;
	long double t = 0.1111111111111111111111111L;
	printf("\nSignificant digits:\n");
	printf("float:       %.25g\n",r);
	printf("double:      %.25g\n",s);
	printf("long double: %.25Lg\n",t);

	return 0;
}
