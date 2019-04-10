#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
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

// Fitting master function
double F(vector* x) {
	double A = vector_get(x,0);
	double T = vector_get(x,1);
	double B = vector_get(x,2);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);
	
	double sum = 0;
	for(int i=0;i<N;i++) sum += pow(A*exp(t[i]/T)+B - y[i],2)/(e[i]*e[i]);
	return sum;
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

void broyden_numerical(
	double f(vector* x),
	vector* x,
	double dx,
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
	
	// B1 - Broyden's update with analytical gradient
	
	fprintf(stderr,"\nBroyden's update with analytical gradient:\n");

	printf("Minimization using Broyden's update, analytical gradient.\n\n");
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("Minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	broyden(f1,set_df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f1(x));
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	broyden(f2,set_df2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f2(x));
	
	// B2 - Broyden's update with numerical gradient
	
	fprintf(stderr,"\nBroyden's update with numerical gradient:\n");
	
	double dx = sqrt(DBL_EPSILON);

	printf("Minimization using Broyden's update, numerical gradient.\n\n");
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("Minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	broyden_numerical(f1,x,dx,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f1(x));
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	broyden_numerical(f2,x,dx,epsilon);
	vector_print(x,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",f2(x));
	
	// B3 - Non-linear least-squares fitting problem
	
	fprintf(stderr,"\nBroyden's update with numerical gradient, fitting problem:\n");
	
	n = 3;
	vector* x_fit = vector_alloc(n);
	
	vector_set(x_fit,0,5);
	vector_set(x_fit,1,-1);
	vector_set(x_fit,2,1);
	printf("Fitting parameters of least squares problem:\n");
	vector_print(x_fit,"x0 =",stdout);
	broyden_numerical(F,x_fit,dx,epsilon);
	vector_print(x_fit,"x_final =",stdout);
	fprintf(stdout,"f(x_final) = %g\n",F(x_fit));
	
	double A = vector_get(x_fit,0);
	double T = vector_get(x_fit,1);
	double B = vector_get(x_fit,2);
	double t0 = 0, tf = 10, dt = (tf-t0)/100;
	FILE* fit_file = fopen("fit_data.txt","w");
	for(double t=t0;t<=tf;t+=dt)
		fprintf(fit_file,"%g %g\n",t,A*exp(t/T)+B);
	
	
	
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
