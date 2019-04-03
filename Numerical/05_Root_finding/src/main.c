#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrix.h"
#include"vector.h"


void f3(vector* x, vector* fx) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(fx,0,2*x1 + 4*x1*(x1*x1+x2-11) + 2*x2*x2 - 14);
	vector_set(fx,1,2*x2 + 4*x2*(x2*x2+x1-7) + 2*x1*x1 - 22);
}

void df3(vector* x, matrix* J) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	matrix_set(J,0,0,12*x1*x1 + 4*x2 - 42);
	matrix_set(J,0,1,4*x1 + 4*x2);
	matrix_set(J,1,0,4*x1 + 4*x2);
	matrix_set(J,1,1,12*x2*x2 + 4*x1 - 26);
}

void f2(vector* x, vector* fx) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(fx,0,2*x1-400*x1*(-(x1*x1)+x2)-2);
	vector_set(fx,1,-200*x1*x1 + 200*x2);
}

void df2(vector* x, matrix* J) {
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	matrix_set(J,0,0,1200*x1*x1 - 400*x2 + 2);
	matrix_set(J,0,1,-400*x1);
	matrix_set(J,1,0,-400*x1);
	matrix_set(J,1,1,200);
}

void f1(vector* x, vector* fx) {
	int A = 10000;
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	vector_set(fx,0,A*x1*x2 - 1);
	vector_set(fx,1,exp(-x1) + exp(-x2) - 1 - 1.0/A);
}

void df1(vector* x, matrix* J) {
	int A = 10000;
	double x1 = vector_get(x,0);
	double x2 = vector_get(x,1);
	matrix_set(J,0,0,A*x2);
	matrix_set(J,0,1,A*x1);
	matrix_set(J,1,0,-exp(-x1));
	matrix_set(J,1,1,-exp(-x2));
}

void newton_with_jacobian(
	void f(vector* x, vector* fx),
	void df(vector* x, matrix* J),
	vector* x,
	double epsilon);

int main() {
	
/*--- A Newton's method with analytic Jacobian and back-tracking linesearch ---*/
	int n = 2;
	vector* x = vector_alloc(n);
	vector* fx = vector_alloc(n);
	double epsilon = 1e-6;
	

	printf("Root finding using analytical jacobian:\n\n");
	vector_set(x,0,8);
	vector_set(x,1,5);
	printf("2 solutions to the first system of equations:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f1,df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f1(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	vector_set(x,0,-2);
	vector_set(x,1,5);
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f1,df1,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f1(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,2);
	vector_set(x,1,1);
	printf("1 minimum of the Rosenbrock valley:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f2,df2,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f2(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_set(x,0,5);
	vector_set(x,1,-5);
	printf("4 minimums of the Himmelblau function:\n");
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	vector_set(x,0,5);
	vector_set(x,1,5);
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	vector_set(x,0,-5);
	vector_set(x,1,5);
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	vector_set(x,0,-5);
	vector_set(x,1,-5);
	vector_print(x,"x0 =",stdout);
	newton_with_jacobian(f3,df3,x,epsilon);
	vector_print(x,"x_final =",stdout);
	f3(x,fx);
	vector_print(fx,"f(x_final) =",stdout);
	
	vector_free(x);
	vector_free(fx);
/*--- B Newton's method with numerical Jacobian and back-tracking linesearch ---*/
/*	
	double dx = 1e-6;
	vector_set(x,0,0);
	vector_set(x,1,0);
	newton(f1,x,dx,epsilon);
	vector_print(x,"Solution to the first system of equations without Jacobian:",stdout);
	vector_set(x,0,0);
	vector_set(x,1,0);
	newton(f2,x,dx,epsilon);
	vector_print(x,"Minimum of the Rosenbrock valley without Jacobian:",stdout);
	vector_set(x,0,0);
	vector_set(x,1,0);
	newton(f3,x,dx,epsilon);
	vector_print(x,"Minimum of the Himmelblau function without Jacobian:",stdout);
*/	
/*--- C Newton's method with refined linesearch ---*/

	return 0;
}
