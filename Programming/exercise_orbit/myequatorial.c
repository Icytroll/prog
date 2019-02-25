#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_equatorial(double t, const double y[], double dydt[], void *params)
{
	double eps = *(double *)params;
	dydt[0] = y[1];
	dydt[1] = 1-y[0]+eps*y[0]*y[0];
	return GSL_SUCCESS;
}

double myequatorial(double t, double EPS, double y0, double yp0)
{
	gsl_odeiv2_system sys;
	sys.function = ode_equatorial;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &EPS;

	gsl_odeiv2_driver *driver;
	double hstart = 0.05, abs = 1e-5, eps = 1e-5;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,
					       gsl_odeiv2_step_rkf45,
					       hstart, abs, eps);

	double t0 = 0;
	double y[] = { y0, yp0 };
	gsl_odeiv2_driver_apply(driver, &t0, t, y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}
