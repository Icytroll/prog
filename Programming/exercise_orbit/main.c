#include<stdio.h>
#include<math.h>
double mylogistic(double x);
double myequatorial(double x, double eps, double y0, double yp0);

int main()
{
	
	/* Integrate the ODE of the logistic function */

	for(double x=0;x<=3;x+=0.05)
		fprintf(stdout,"%g %g %g\n",x,mylogistic(x),1/(1+exp(-x)));
	

	/* Integrate the ODE of equatorial motion of a planet */

	double eps[3] = {0,0,0.1};
	double y0[3]  = {1,1,1,};
	double yp0[3] = {0,0.5,0.5};

	for(double phi=0;phi<=10*2*M_PI;phi+=0.05)
		fprintf(stderr,"%g %g %g %g\n",phi,myequatorial(phi,eps[0],y0[0],yp0[0])
										  ,myequatorial(phi,eps[1],y0[1],yp0[1])
										  ,myequatorial(phi,eps[2],y0[2],yp0[2]));

	return 0;
}
