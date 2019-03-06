#include<stdio.h>
#include<stdlib.h>

void rosenbrock(double x0, double y0);
void LS_fit(int n, double* t, double* y, double* e);

int main(int argc, char** argv) {
	
	/* Find minimum of the Rosenbrock function using GSL minimization functions */
	
	double x0, y0;
	if(argc<2) {
		x0 = -3;
		y0 = -3;
	}
	else {
		x0 = atof(argv[1]);
		y0 = atof(argv[2]);
	}	

	rosenbrock(x0,y0);
	
	/* Least squares fit of experimental data */
	
	int n = atoi(argv[3]);
	double t[n], y[n], e[n];
	for(int i=0;i<n;i++)
		scanf("%lg %lg %lg",t+i,y+i,e+i);	
		
	LS_fit(n,t,y,e);
	
	return 0;
}
