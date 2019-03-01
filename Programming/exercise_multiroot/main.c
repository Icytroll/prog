#include<stdio.h>

int rosen_extremum(double x0, double y0);

int main() {

	int flag;
	flag = rosen_extremum(3.1,-1);	
	
	fprintf(stderr,"flag = %d\n",flag);

	return 0;
}
