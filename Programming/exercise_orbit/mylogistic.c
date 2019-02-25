#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_logistic(double t, const double y[], double dydt[], void *params)
{
	dydt[0] = y[0]*(1-y[0]);
	return GSL_SUCCESS;
}

double mylogistic(double t)
{
	gsl_odeiv2_system sys;
	sys.function = ode_logistic;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = NULL;

	gsl_odeiv2_driver *driver;
	double hstart = 0.05, abs = 1e-5, eps = 1e-5;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,
					       gsl_odeiv2_step_rkf45,
					       hstart, abs, eps);

	double t0 = 0;
	double y[] = { 0.5 };
	gsl_odeiv2_driver_apply(driver, &t0, t, y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}
