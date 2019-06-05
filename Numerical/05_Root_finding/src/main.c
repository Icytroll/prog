#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"matrix.h"
#include"vector.h"

/*--- FUNCTION DECLARATIONS ---*/

// ROOT FUNCTIONS AND THEIR DERIVATIVES

void f1  (vector* x, vector* fx); // System of equations
void df1 (vector* x, matrix* J);
void f2  (vector* x, vector* fx); // Rosenbrock Valley
void df2 (vector* x, matrix* J);
void f3  (vector* x, vector* fx); // Himmelblau Function
void df3 (vector* x, matrix* J);

void f1_gsl   (const gsl_vector* x,void* params,gsl_vector* f);
void df1_gsl  (const gsl_vector* x,void* params,gsl_matrix* J);
void fdf1_gsl (const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J);
void f2_gsl   (const gsl_vector* x,void* params,gsl_vector* f);
void df2_gsl  (const gsl_vector* x,void* params,gsl_matrix* J);
void fdf2_gsl (const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J);
void f3_gsl   (const gsl_vector* x,void* params,gsl_vector* f);
void df3_gsl  (const gsl_vector* x,void* params,gsl_matrix* J);
void fdf3_gsl (const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J);

// ROOT SOLVERS

void newton_with_jacobian(
	void f(vector* x, vector* fx),
	void df(vector* x, matrix* J),
	vector* x,
	double epsilon);

void newton(
	void f(vector* x, vector* fx),
	vector* x,
	double dx,
	double epsilon);

void newton_quadline(
	void f(vector* x, vector* fx),
	void df(vector* x, matrix* J),
	vector* x,
	double epsilon);

void gsl_rootfinder(
	void f(const gsl_vector* x, void* params, gsl_vector* f),
    void df(const gsl_vector* x, void* params, gsl_matrix* J),
    void fdf(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J),
    gsl_vector* x,
    double epsilon);

/*--- MAIN PROGRAM ---*/

int main() {
	
	fprintf(stderr,"ITERATIONS, FUNCTION CALLS:\n\n");
	fprintf(stderr,"(System of equations)\n");
	fprintf(stderr,"(Rosenbrock Valley)\n");
	fprintf(stderr,"(Himmelblau function)\n");
	
	// A - Newton's method with analytic Jacobian and back-tracking linesearch
	
	fprintf(stderr,"\nNewton with analytic jacobian:\n");
	
	int n = 2;
	vector* x = vector_alloc(n);
	vector* fx = vector_alloc(n);
	double epsilon = 1e-10;
	
	printf("Root finding using analytical jacobian:\n\n");
	vector_set(x,0,3);
	vector_set(x,1,6);
	printf("1 of the solutions to the first system of equations:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f1,df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f1(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("Minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f2,df2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f2(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	// B - Newton's method with numerical Jacobian and back-tracking linesearch
	// Also testing GSL root multiroot finder
	
	fprintf(stderr,"\nNewton with numerical jacobian:\n");
	
	double dx = sqrt(DBL_EPSILON);
	gsl_vector* x_gsl = gsl_vector_alloc(n);

	printf("Root finding using numerical jacobian\n\n");
	
	vector_set(x,0,3);
	vector_set(x,1,6);
	printf("1 of the solutions to the first system of equations:\n");
	vector_print(x,"x0 =",stdout);
	newton(f1,x,dx,epsilon);
	vector_print(x,"x_final =",stdout);
	f1(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	vector_print(x,"x0 =",stdout);
	printf("Minimum of the Rosenbrock valley:\n");
	newton(f2,x,dx,epsilon);
	vector_print(x,"x_final =",stdout);
	f2(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton(f3,x,dx,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	fprintf(stderr,"\nGSL rootfinder (only iterations):\n");

	gsl_vector_set(x_gsl,0,3);
	gsl_vector_set(x_gsl,1,6);
	gsl_rootfinder(f1_gsl,df1_gsl,fdf1_gsl,x_gsl,epsilon);
	
	gsl_vector_set(x_gsl,0,2);
	gsl_vector_set(x_gsl,1,1);
	gsl_rootfinder(f2_gsl,df2_gsl,fdf2_gsl,x_gsl,epsilon);
	
	gsl_vector_set(x_gsl,0,2.5);
	gsl_vector_set(x_gsl,1,1.5);
	gsl_rootfinder(f3_gsl,df3_gsl,fdf3_gsl,x_gsl,epsilon);
	

	// C - Newton's method with refined linesearch
	
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
	
	vector_set(x,0,2.5);
	vector_set(x,1,1.5);
	printf("1 of the minimum of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton_quadline(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	// Clean up

	vector_free(x);
	vector_free(fx);
	
	return 0;
}
