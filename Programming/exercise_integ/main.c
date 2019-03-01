#include<stdio.h>
#include<math.h>
double logSqrt(double a, double b);
double E_numerical(double alpha);

double E_analytical(double alpha) {
	return ((sqrt(M_PI)*(alpha*alpha+1))/(4*pow(alpha,1.5)))/(sqrt(M_PI)/sqrt(alpha));
}

int main() {

	/* Numerically integrate log(x)/sqrt(x) from 0<=x<=1 */

	double result = logSqrt(0,1);
	fprintf(stdout,"Integration result of log(x)/sqrt(x): %g\n",result);	
	
	/* Numerically integrate the expectation value E(alpha) of a Hamiltonian H */
	
	for(double alpha=0.001;alpha<=5;alpha+=0.001)
		fprintf(stderr,"%g %g %g\n",alpha,E_numerical(alpha),E_analytical(alpha));	
	
	return 0;
}
