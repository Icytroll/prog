#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include"qspline.h"
#include"cspline.h"

double linterp(int n, double *x, double *y, double z);
double linterp_integ(int n, double *x, double *y, double z);

int main() {
	
	// Generate tabulated data from a sine function
	FILE * tabdata = fopen("tabdata.txt","w");
	
	int n = 5;
	double xlim = 2*M_PI,dx = xlim/(n-1);
	double x[n],y[n];
	for(int i=0;i<n;i++) {
		x[i] = dx*i;
		y[i] = sin(dx*i);
		fprintf(tabdata,"%g %g\n",x[i],y[i]);
	}	
	
	fclose(tabdata);
	
//------------------------------------------------------     
	FILE * fid = fopen("data.txt","w");
	
	/* A - Linear splines */
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g %g %g %g\n",z,linterp(n,x,y,z),sin(z)
									  ,linterp_integ(n,x,y,z),1-cos(z));
	
	fprintf(fid,"\n\n");

	/* B - Quadratic splines */
	qspline *s = qspline_alloc(n,x,y);
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g %g %g %g %g %g\n",z,qspline_eval(s,z),sin(z)
											,qspline_deriv(s,z),cos(z)
											,qspline_integ(s,z),1-cos(z));
	qspline_free(s);
	
	fprintf(fid,"\n\n");
	
	/* C - Cubic splines */
	// GSL cubic spline
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,n);
	gsl_spline_init(spline,x,y,n);
	
	// own formulation
	cspline *S = cspline_alloc(n,x,y);
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g %g %g %g %g %g %g\n",z,cspline_eval(S,z),sin(z)
											,gsl_spline_eval(spline,z,acc)
											,cspline_deriv(S,z),cos(z)
											,cspline_integ(S,z),1-cos(z));
	cspline_free(S);
	
	fclose(fid);

	return 0;
}
