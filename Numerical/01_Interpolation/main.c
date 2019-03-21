#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"qspline.h"

double linterp(int n, double *x, double *y, double z);
double linterp_integ(int n, double *x, double *y, double z);



int main(int argc, char** argv) {
	
	// Generate tabulated data from a sine function
	int n = 5;
	double xlim = 2*M_PI,dx = xlim/(n-1);
	double x[n],y[n];
	for(int i=0;i<n;i++) {
		x[i] = dx*i;
		y[i] = sin(dx*i);
		printf("x=%g y=%g\n",x[i],y[i]);
	}	
	
	FILE * fid = fopen("data.txt","w");
	
	/* A - Linear splines */
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g %g %g %g\n",z,linterp(n,x,y,z),sin(z)
									  ,linterp_integ(n,x,y,z),1-cos(z));
	
	fprintf(fid,"\n\n");

	/* B - Quadratic splines */
	printf("Se mig 1\n");
	qspline *s = qspline_alloc(n,x,y);
	printf("Se mig 2\n");
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g\n",z,qspline_eval(s,z));
		/*fprintf(fid,"%g %g %g %g %g %g %g\n",z,qspline_eval(s,z),sin(z)
											,qspline_deriv(s,z),cos(z)
											,qspline_integ(s,z),1-cos(z))*/;
	printf("Se mig 3\n");
	fprintf(fid,"\n\n");
	
	/* C - Cubic splines 
	for(double z=0;z<=2*M_PI;z+=M_PI/50)
		fprintf(fid,"%g %g %g %g %g %g %g\n",z,cspline_eval(s,z),sin(z)
											,cspline_deriv(s,z),cos(z)
											,cspline_integ(s,z),-cos(z));
	*/
	return 0;
}
