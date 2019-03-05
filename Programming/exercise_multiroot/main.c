#include<stdio.h>

int rosen_extremum(double x0, double y0);
double s_wave_root(double eps0, double rmax);
double fe(double eps, double rmax, int saveData);

int main() {
	
	/* Find the extremum of the Rosenbrock function */

	int flag;
	flag = rosen_extremum(3.1,-1);	
	
	fprintf(stderr,"Rosenbrock flag = %d\n",flag);

	/* Use the Shooting Method to find the lowest root of the equation M(eps) = 0 */
	
	double eps = s_wave_root(-1,8);
	fprintf(stderr,"Lowest root of M(eps) = %g\n",eps);
	
	for(int r=2;r<=10;r+=2) fe(eps,r,1);
		
	return 0;
}
