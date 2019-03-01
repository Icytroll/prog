#include<stdio.h>
#include<math.h>
double logSqrt(double a, double b);
double expecValue(double alpha);

int main() {

	/* Numerically integrate log(x)/sqrt(x) from 0<=x<=1 */

	double result = logSqrt(0,1);
	fprintf(stdout,"Integration result of log(x)/sqrt(x): %g\n",result);	
	
	/* Numerically integrate the expectation value E(alpha) of a Hamiltonian H */
	
	for(double alpha=0.001;alpha<=5;alpha+=0.001)
		fprintf(stderr,"%g %g\n",alpha,expecValue(alpha));	
	
	return 0;
}
