#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double gsl_integ(double f(double x, void* params), double a, double b, double acc, double eps) {
	int limit = 10000;
	gsl_integration_workspace* ws = gsl_integration_workspace_alloc(limit);
	gsl_function F;
	F.function = f;
	F.params = NULL;
	double result, error;
	int flag = gsl_integration_qags(&F,a,b,acc,eps,limit,ws,&result,&error);
	gsl_integration_workspace_free(ws);
	if(flag!=GSL_SUCCESS) return NAN;
	return result;
}
