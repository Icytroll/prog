#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double numerator (double x, void * params) {
	double alpha = *(double*)params;
	double f = (-(alpha*alpha)*(x*x)/2 + alpha/2 + (x*x)/2)*exp(-alpha*(x*x));
	return f;
}

double denominator (double x, void * params) {
	double alpha = *(double*)params;
	double f = exp(-alpha*(x*x));
	return f;
}


double E_numerical(double alpha){
	
	int limit=1000;
	gsl_integration_workspace * wN;
	wN = gsl_integration_workspace_alloc (limit);

	gsl_integration_workspace * wD;
	wD = gsl_integration_workspace_alloc (limit);

	double params = alpha;

	gsl_function N;
	N.function = numerator;
	N.params = (void*)&params;

	gsl_function D;
	D.function = denominator;
	D.params = (void*)&params;
	
	double result_N, error_N, epsabs = 1e-3, epsrel = 1e-3;
	int flag_N = gsl_integration_qagi
		(&N, epsabs, epsrel, limit, wN, &result_N, &error_N);
	
	double result_D, error_D;
	int flag_D = gsl_integration_qagi
		(&D, epsabs, epsrel, limit, wD, &result_D, &error_D);
	
	gsl_integration_workspace_free(wN);
	gsl_integration_workspace_free(wD);

	if(flag_N!=GSL_SUCCESS || flag_D!=GSL_SUCCESS) return NAN;
	
	return result_N/result_D;
}
