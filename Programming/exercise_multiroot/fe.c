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

double fe(double eps, double rmax, int saveData)
{
	double rmin = 1e-3;
	if(rmax<rmin) return rmax-rmax*rmax;

	gsl_odeiv2_system sys;
	sys.function = s_wave_ode;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)&eps;

	gsl_odeiv2_driver *driver;

	double r = rmin;
	double dr = (rmax-rmin)/100;
	double ABS = 1e-12, EPS = 1e-12;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, dr, ABS, EPS);

	double y[] = {rmin-rmin*rmin,1-2*rmin};
	
	FILE * stream;
	stream = fopen("sData.txt","a");

	if(saveData) fprintf(stream,"# Rmax = %g\n",rmax);

	for(double ri=r+dr;ri<=rmax;ri+=dr) {
		gsl_odeiv2_driver_apply(driver, &r, ri, y);
		if(saveData) fprintf(stream,"%g %g %g\n",ri,y[0],ri*exp(-ri));
	}
	
	if(saveData) fprintf(stream,"\n");
	fclose(stream);

	gsl_odeiv2_driver_free(driver);

	return y[0];
}
