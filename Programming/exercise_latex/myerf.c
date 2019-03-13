#include<math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_efunc(double t, const double y[], double dydt[], void *params)
{
	dydt[0] = 2/sqrt(M_PI)*exp(-(t*t));
	return GSL_SUCCESS;
}

double myerf(double x, double dx) {
	
	// define integration parameters
	if(x<0) 
		return -myerf(-x,dx);
	
	// define ode system
	gsl_odeiv2_system sys;
	sys.function = ode_efunc;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = NULL;
	
	// define ode driver
	gsl_odeiv2_driver *driver;
	double abs = 1e-8, eps = 1e-8;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,
					       gsl_odeiv2_step_rkf45,
					       dx, abs, eps);
	
	// apply the driver
	double t0 = 0;
	double y[] = { 0 };
	gsl_odeiv2_driver_apply(driver, &t0, x, y);
	
	// cleanup
	gsl_odeiv2_driver_free(driver);
	
	return y[0];
}
