#include<stdio.h>

int rosen_extremum(double x0, double y0);
double s_wave_root(double eps0, double rmax);
double fe(double eps, double rmax, FILE* stream);

int main() {
	
	/* Find the extremum of the Rosenbrock function */

	int flag;
	flag = rosen_extremum(3.1,-1);	
	
	fprintf(stderr,"Rosenbrock flag = %d\n",flag);

	/* Use the Shooting Method to find the lowest root of the equation M(eps) = 0 */
	

	// investigate different values of Rmax
	double rmax,eps;
	FILE* stream=fopen("sData.txt","w");

	rmax=2;
	eps = s_wave_root(-1,rmax);
	fprintf(stderr,"Lowest root of M(eps) = %g\n",eps);
	fe(eps,rmax,stream);

	rmax=4;
	eps = s_wave_root(-1,rmax);
	fprintf(stderr,"Lowest root of M(eps) = %g\n",eps);
	fe(eps,rmax,stream);

	rmax=6;
	eps = s_wave_root(-1,rmax);
	fprintf(stderr,"Lowest root of M(eps) = %g\n",eps);
	fe(eps,rmax,stream);

	rmax=8;
	eps = s_wave_root(-1,rmax);
	fprintf(stderr,"Lowest root of M(eps) = %g\n",eps);
	fe(eps,rmax,stream);

	return 0;
}
