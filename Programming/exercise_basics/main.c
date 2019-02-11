#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<float.h>
#include<stdlib.h>

int equal(double a, double b, double tau, double epsilon);

int main(int argc, char** argv) {
	
	// Exit program if not enough inputs for "equal" function
	if (argc < 5) {
		printf("Please provide 4 inputs for the equal function!\n");
		exit(0);
	}

	int i = 1;

	// Find max representable integer
	
	while (i+1>i) {i++;}
	printf("Max integer, while:    %i\n",i);

	for (i=1;i+1>i;i++) {
		if (i+2<=i) printf("Max integer, for:      %i\n",i+1);
	}

	i = 1;
	do (i++); while (i+1>i);
	printf("Max integer, do while: %i\n",i);

	printf("INT_MAX =              %i\n\n",INT_MAX);
	
	// Find min representable integer

	i = 1;
	while (i-1<i) {i--;}
	printf("Min integer, while:    %i\n",i);

	for (i=1;i-1<i;i--) {
		if (i-2>=i) printf("Min integer, for:      %i\n",i-1);
	}

	i = 1;
	do (i--); while (i-1<i);
	printf("Max integer, do while: %i\n",i);

	printf("INT_MIN =              %i\n\n",INT_MIN);
	
	// Find machine epsilon
	// Float
	
	float x = 1;
	double y = 1;
	long double z = 1;

	while (1+x!=1) {x/=2;}
	x*=2;
	printf("Machine epsilon, float, while:    %g\n",x);

	for (x=1;1+x!=1;x/=2) {
		if (1+x/2==1) printf("Machine epsilon, float, for:      %g\n",x);
	}

	x = 1;
	do (x/=2); while (1+x!=1);
	x*=2;
	printf("Machine epsilon, float, do while: %g\n",x);
	
	printf("FLT_EPSILON =                     %g\n\n",FLT_EPSILON);

	// Double

	while (1+y!=1) {y/=2;}
	y*=2;
	printf("Machine epsilon, double, while:    %g\n",y);

	for (y=1;1+y!=1;y/=2) {
		if (1+y/2==1) printf("Machine epsilon, double, for:      %g\n",y);
	}

	y = 1;
	do (y/=2); while (1+y!=1);
	y*=2;
	printf("Machine epsilon, double, do while: %g\n",y);
	
	printf("DBL_EPSILON =                      %g\n\n",DBL_EPSILON);
	
	// Long double

	while (1+z!=1) {z/=2;}
	z*=2;
	printf("Machine epsilon, float, while:    %Lg\n",z);

	for (z=1;1+z!=1;z/=2) {
		if (1+z/2==1) printf("Machine epsilon, float, for:      %Lg\n",z);
	}

	z = 1;
	do (z/=2); while (1+z!=1);
	z*=2;
	printf("Machine epsilon, float, do while: %Lg\n",z);
	
	printf("LDBL_EPSILON =                    %Lg\n\n",LDBL_EPSILON);
	
	// Float and double sums
	
	int max = INT_MAX/2;
	float  sum_up_float  = 0, sum_down_float  = 0;
	double sum_up_double = 0, sum_down_double = 0;

	for (int i=1;i<=max;i++) {
		sum_up_float    += 1.0f/i;
		sum_down_float  += 1.0f/(max-i+1);
		sum_up_double   += 1.0f/i;
		sum_down_double += 1.0f/(max-i+1);
	}
	
	printf("sum_up_float    = %g\n",sum_up_float);
	printf("sum_down_float  = %g\n",sum_down_float);
	printf("\nThe difference between the above two sums is that sum_up_float fills up its variable with lots of numbers before it gets to the small stuff, at which point the precision of the float is all used up.\nThe float precision is not as much a problem for sum_down, as it includes the small stuff first and therefore keeps that precision until the end of the sum.\n\n");

	printf("The float sums do seem to converge to the above numbers, as we can't add anything extra to the sum once max becomes large enough.\n\n");
	printf("sum_up_double   = %g\n",sum_up_double);
	printf("sum_down_double = %g\n",sum_down_double);
	printf("The double sums are equal since there in this case is plenty of precision to go around, as INT_MAX only uses ~10 digits of precision and a double has precision up to 16 digits.\n\n");

	// Equal function
	
	// extract input
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double tau = atof(argv[3]);
	double epsilon = atof(argv[4]);

	int res = equal(a,b,tau,epsilon);
	printf("a and b are equal: %i\n",res);	
	
	return 0;
}


































