#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf.h>

double myerf(double x, double dx);

int main(int argc, char** argv) {
	
	double a,b,dx;
	if(argc<2) a = -3;
	else a = atof(argv[1]);
	if(argc<3) b = 3;
	else b = atof(argv[2]);
	if(argc<4) dx = 0.1;
	else dx = atof(argv[3]);

	for(double i=a;i<=b;i+=dx)
		printf("%g %g %g\n",i,myerf(i,dx),gsl_sf_erf(i));
	
	return 0;
}
