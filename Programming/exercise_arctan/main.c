#include<stdio.h>
#include<math.h>

double myarctan(double x);

int main() {
	
	for(double x=-20;x<=20;x+=0.01)
		printf("%g %g %g\n",x,myarctan(x),atan(x));

	return 0;
}
