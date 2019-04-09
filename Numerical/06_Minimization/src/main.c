#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"matrix.h"
#include"vector.h"

/*--- FUNCTION DECLARATIONS ---*/

// TEST-FUNCTIONS AND THEIR DERIVATIVES/HESSIAN

// Rosenbrocks Valley
double f1 (vector* x) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	return (1-x1)*(1-x1) + 100*(x2-x1*x1)*(x2-x1*x1);
}

void set_df1 (vector* x, vector* df) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(df,0,2*x1 - 400*x1*(-(x1*x1)+x2) - 2);
	vector_set(df,1,-200*x1*x1 + 200*x2);
}

void set_dfH1 (vector* x, vector* df, matrix* H) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(df,0,2*x1 - 400*x1*(-(x1*x1)+x2) - 2);
	vector_set(df,1,-200*x1*x1 + 200*x2);
	matrix_set(H,0,0,1200*x1*x1 - 400*x2 + 2);
	matrix_set(H,0,1,-400*x1);
	matrix_set(H,1,0,-400*x1);
	matrix_set(H,1,1,200);
}

// Himmelblau Function
double f2 (vector* x) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	return pow(x1*x1+x2-11,2) + pow(x1+x2*x2-7,2);
}

void set_df2 (vector* x, vector* df) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(df,0,2*x1 + 4*x1*(x1*x1+x2-11) + 2*x2*x2 - 14);
	vector_set(df,1,2*x2 + 4*x2*(x2*x2+x1-7) + 2*x1*x1 - 22);
}

void set_dfH2 (vector* x, vector* df, matrix* H) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(df,0,2*x1 + 4*x1*(x1*x1+x2-11) + 2*x2*x2 - 14);
	vector_set(df,1,2*x2 + 4*x2*(x2*x2+x1-7) + 2*x1*x1 - 22);
	matrix_set(H,0,0,12*x1*x1 + 4*x2 - 42);
	matrix_set(H,0,1,4*x1 + 4*x2);
	matrix_set(H,1,0,4*x1 + 4*x2);
	matrix_set(H,1,1,12*x2*x2 + 4*x1 - 26);
}

// MINIMIZATION FUNCTIONS

void newton(
	double f(vector* x),
	void set_dfH(vector* x, vector* df, matrix* H),
	vector* x,
	double epsilon);

void broyden(
	double f(vector* x),
	void set_df(vector* x, vector* df),
	vector* x,
	double epsilon);

/*--- MAIN PROGRAM ---*/

int main() {
	
	fprintf(stderr,"ITERATIONS:\n\n");
	fprintf(stderr,"(Rosenbrock Valley)\n");
	fprintf(stderr,"(Himmelblau function)\n");
	
	// A - Newton's method with analytic gradient and hessian, and back-tracking linesearch
	
	fprintf(stderr,"\nNewton with analytic gradient and hessian:\n");
	
	int n = 2;
	vector* x = vector_alloc(n);
	double epsilon = 1e-10;
	
	printf("Minimization using analytical gradient and hessian:\n\n");
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("Minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	newton(f1,set_dfH1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f1(x));
		
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton(f2,set_dfH2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f2(x));
	
	// B - Broyden's update with analytical gradient
	
	fprintf(stderr,"\nBroyden's update with analytical gradient:\n");
	
	double dx = sqrt(DBL_EPSILON);

	printf("Minimization using Broyden's update.\n\n");
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("Minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	broyden(f1,set_df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f1(x));
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	broyden(f2,set_df2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f2(x));
	
	// C - Newton's method with refined linesearch
	/*
	fprintf(stderr,"\nNewton with refined linesearch:\n");
	
	printf("Root finding using refined linesearch\n\n");
	
	vector_set(x,0,3);
	vector_set(x,1,6);
	printf("1 of the solutions to the first system of equations:\n");
	vector_print(x,"x0 =",stdout);
	newton_quadline(f1,df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f1(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	vector_print(x,"x0 =",stdout);
	printf("Minimum of the Rosenbrock valley:\n");
	newton_quadline(f2,df2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f2(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,5);
	vector_set(x,1,5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton_quadline(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	*/
	// Clean up

	vector_free(x);
	
	return 0;
}
