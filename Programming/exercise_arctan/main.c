#include<stdio.h>
#include<math.h>

double myarctan(double x);

int main() {
	
	for(double x=-3;x<=3;x+=0.1)
		printf("%g %g %g\n",x,myarctan(x),atan(x));

	return 0;
}
