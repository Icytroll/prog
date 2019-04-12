#include"vector.h"
#include"matrix.h"


void driver(
    double* t,
    double b,
    double* h,
    vector* yt,
    double acc, double eps,
    void stepper(
        double t, double h, vector* yt,
        void f(double t, vector* y, vector* dydt, void* params),
        vector* yth, vector* err),
    void f(double t, vector* y, vector* dydt, void* params)) {
	
	
}
