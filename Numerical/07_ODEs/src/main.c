#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include"matrix.h"
#include"vector.h"

/*--- TEST FUNCTIONS ---*/

// Equation of motion
void EoM(double t, vector* y, vector* dydt, void* params) {
	double y1 = vector_get(y,0);
	double y2 = vector_get(y,1);
	double* p = (double*)params;
	double m = p[0], c = p[1], k = p[2];
	vector_set(dydt,0,y2);
	vector_set(dydt,1,-c/m*y2 - k/m*y1);
}

// Van der Pol oscillator
void VanPol(double t, vector* y, vector* dydt, void* params) {
	double y1 = vector_get(y,0);
	double y2 = vector_get(y,1);
	double mu = *(double*)params;
	vector_set(dydt,0,y2);
	vector_set(dydt,1,mu*(1-y1*y1)*y2 - y1);
}

/*--- FUNCTION DECLARATIONS ---*/

// Runge-Kutta 45 stepper
void rkstep45(
	double t,
	double h,
	vector* yt,
	void f(double t, vector* y, vector* dydt, void* params),
	void* params,
	vector* yth,
	vector* err);

// ODE driver
int driver(
	double* t,
	double b,
	double* h,
	vector* yt,
	double acc, double eps,
	void stepper(
		double t, double h, vector* yt,
		void f(double t, vector* y, vector* dydt, void* params),
		void* params, vector* yth, vector* err),
	void f(double t, vector* y, vector* dydt, void* params),
	void* params,
	matrix* data);

/*--- MAIN PROGRAM ---*/

int main() {
	
	// A - Embedded Runge-Kutta ODE integrator
	// b - Storing the path

	int n = 2, m = 10000;
	vector* yt = vector_alloc(n);
	matrix* data = matrix_alloc(m,n+1);
	double t0 = 0, h0 = 0.01;
	double* t = &t0;
	double* h = &h0;
	double b = 20, acc = 1e-7, eps = 1e-7, params[3] = {1.0,0.1,1.0};
	vector_set(yt,0,1);
	vector_set(yt,1,0);
	
	printf("Mass-damper-spring, m = 1, c = 0.1, k = 1:\n");
	vector_print(yt,"y0 =",stdout);
	int N = driver(t,b,h,yt,acc,eps,rkstep45,EoM,(void*)params,data);
	vector_print(yt,"yf =",stdout);
	
	for(int i=0;i<N;i++) {
		for(int j=0;j<n+1;j++) {
			fprintf(stderr,"%g ",matrix_get(data,i,j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n\n");

	vector_set(yt,0,1);
	vector_set(yt,1,0);
	params[1] = 1;
	*t = 0; *h = 0.01;
	printf("Mass-damper-spring, m = 1, c = 1, k = 1:\n");
	vector_print(yt,"y0 =",stdout);
	N = driver(t,b,h,yt,acc,eps,rkstep45,EoM,(void*)params,data);
	vector_print(yt,"yf =",stdout);
	for(int i=0;i<N;i++) {
		for(int j=0;j<n+1;j++) {
			fprintf(stderr,"%g ",matrix_get(data,i,j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n\n");
	
	vector_set(yt,0,1);
	vector_set(yt,1,0);
	params[1] = 10;
	*t = 0; *h = 0.01;
	printf("Mass-damper-spring, m = 1, c = 10, k = 1:\n");
	vector_print(yt,"y0 =",stdout);
	N = driver(t,b,h,yt,acc,eps,rkstep45,EoM,(void*)params,data);
	vector_print(yt,"yf =",stdout);
	for(int i=0;i<N;i++) {
		for(int j=0;j<n+1;j++) {
			fprintf(stderr,"%g ",matrix_get(data,i,j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n\n");
	
	


	vector_set(yt,0,2);
	vector_set(yt,1,0);
	*t = 0; *h = 0.01;
	double mu = 1;
	double* Params = &mu;
	printf("Van der Pol equation, mu = 1:\n");
	vector_print(yt,"y0 =",stdout);
	N = driver(t,b,h,yt,acc,eps,rkstep45,VanPol,(void*)Params,data);
	vector_print(yt,"yf =",stdout);
	for(int i=0;i<N;i++) {
		for(int j=0;j<n+1;j++) {
			fprintf(stderr,"%g ",matrix_get(data,i,j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n\n");
	
	
	
	// C - A definite integral as an ODE 
	
	
	matrix_free(data);
	vector_free(yt);
	return 0;
}
