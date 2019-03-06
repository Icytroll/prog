#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>

int s_wave_ode(double r, const double y[], double dydr[], void *params)
{
	double eps = *(double*)params;
	dydr[0] = y[1];
	dydr[1] = 2*y[0]*(-1/r-eps);
	return GSL_SUCCESS;
}

double fe(double eps, double rmax, FILE* stream)
{
	double rmin = 1e-3;
	if(rmax<rmin) return rmax-rmax*rmax;
	
	// define ode system
	gsl_odeiv2_system sys;
	sys.function = s_wave_ode;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)&eps;

	gsl_odeiv2_driver *driver;

	double r = rmin;
	double dr = (rmax-rmin)/100;
	double ABSEPS = 1e-12, RELEPS = 1e-12;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, dr, ABSEPS, RELEPS);

	double y[] = {rmin-rmin*rmin,1-2*rmin};
	
	// prepare to write data
	
	if(stream!=NULL) fprintf(stream,"# Rmax = %g\n",rmax);
	
	// manually step the solution 
	if(stream!=NULL){
		for(double ri=r+dr;ri<=rmax;ri+=dr) {
			gsl_odeiv2_driver_apply(driver, &r, ri, y);
			fprintf(stream,"%g %g %g\n",r,y[0],r*exp(-r));
		}
		if(stream!=NULL) fprintf(stream,"\n\n");
	}

	r=rmin;
	gsl_odeiv2_driver_apply(driver, &r, rmax, y);
	
	
	// cleanup
	gsl_odeiv2_driver_free(driver);

	return y[0];
}
