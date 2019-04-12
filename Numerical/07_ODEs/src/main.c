#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include"matrix.h"
#include"vector.h"

/*--- FUNCTION DECLARATIONS ---*/

// Equation of motion
void EoM(double t, vector* y, vector* dydt, void* params) {
	double y1 = vector_get(y,0);
	double y2 = vector_get(y,1);
	double m = *(double*)params[0];
	double c = *(double*)params[1];
	double k = *(double*)params[2];
	vector_set(dydt,0,y2);
	vector_set(dydt,1,-c/m*y2 - k/m*y1);
}

// Runge-Kutta 45 stepper
void rkstep45(
	double t,
	double h,
	vector* yt,
	void f(double t, vector* y, vector* dydt),
	vector* yth,
	vector* err
)

// ODE driver
void driver(
	double* t,
	double b,
	double* h,
	vector* yt,
	double acc, double eps,
	void stepper(
		double t, double h, vector* yt,
		void f(double t, vector* y, vector* dydt),
		vector* yth, vector* err),
	void f(double t, vector* y, vector* dydt)
)

/*--- MAIN PROGRAM ---*/

int main() {
	
	// A - Embedded Runge-Kutta ODE integrator
	
	
	
	// B - Storing the path
	
	// C - A definite integral as an ODE 

	return 0;
}
