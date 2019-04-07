#include<math.h>
#include"vector.h"
#include"matrix.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_errno.h>

// ROOT FUNCTIONS AND THEIR DERIVATIVES

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

int f1_gsl(const gsl_vector* x,void* params,gsl_vector* f) {
	int A = 10000;
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,A*x1*x2 - 1);
    gsl_vector_set(f,1,exp(-x1) + exp(-x2) - 1 - 1.0/A);
	return GSL_SUCCESS;
}

int df1_gsl(const gsl_vector* x,void* params,gsl_matrix* J) {
    int A = 10000;
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_matrix_set(J,0,0,A*x2);
    gsl_matrix_set(J,0,1,A*x1);
    gsl_matrix_set(J,1,0,-exp(-x1));
    gsl_matrix_set(J,1,1,-exp(-x2));
	return GSL_SUCCESS;
}

int fdf1_gsl(const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J) {
	int A = 10000;
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,A*x1*x2 - 1);
    gsl_vector_set(f,1,exp(-x1) + exp(-x2) - 1 - 1.0/A);
    gsl_matrix_set(J,0,0,A*x2);
    gsl_matrix_set(J,0,1,A*x1);
    gsl_matrix_set(J,1,0,-exp(-x1));
    gsl_matrix_set(J,1,1,-exp(-x2));
	return GSL_SUCCESS;
}

int f2_gsl(const gsl_vector* x,void* params,gsl_vector* f) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,2*x1-400*x1*(-(x1*x1)+x2)-2);
    gsl_vector_set(f,1,-200*x1*x1 + 200*x2);
	return GSL_SUCCESS;
}

int df2_gsl(const gsl_vector* x,void* params,gsl_matrix* J) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_matrix_set(J,0,0,1200*x1*x1 - 400*x2 + 2);
    gsl_matrix_set(J,0,1,-400*x1);
    gsl_matrix_set(J,1,0,-400*x1);
    gsl_matrix_set(J,1,1,200);
	return GSL_SUCCESS;
}

int fdf2_gsl(const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,2*x1-400*x1*(-(x1*x1)+x2)-2);
    gsl_vector_set(f,1,-200*x1*x1 + 200*x2);
    gsl_matrix_set(J,0,0,1200*x1*x1 - 400*x2 + 2);
    gsl_matrix_set(J,0,1,-400*x1);
    gsl_matrix_set(J,1,0,-400*x1);
    gsl_matrix_set(J,1,1,200);
	return GSL_SUCCESS;
}

int f3_gsl(const gsl_vector* x,void* params,gsl_vector* f) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,2*x1 + 4*x1*(x1*x1+x2-11) + 2*x2*x2 - 14);
    gsl_vector_set(f,1,2*x2 + 4*x2*(x2*x2+x1-7) + 2*x1*x1 - 22);
	return GSL_SUCCESS;
}

int df3_gsl(const gsl_vector* x,void* params,gsl_matrix* J) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_matrix_set(J,0,0,12*x1*x1 + 4*x2 - 42);
    gsl_matrix_set(J,0,1,4*x1 + 4*x2);
    gsl_matrix_set(J,1,0,4*x1 + 4*x2);
    gsl_matrix_set(J,1,1,12*x2*x2 + 4*x1 - 26);
	return GSL_SUCCESS;
}

int fdf3_gsl(const gsl_vector* x,void* params,gsl_vector* f,gsl_matrix* J) {
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    gsl_vector_set(f,0,2*x1 + 4*x1*(x1*x1+x2-11) + 2*x2*x2 - 14);
    gsl_vector_set(f,1,2*x2 + 4*x2*(x2*x2+x1-7) + 2*x1*x1 - 22);
    gsl_matrix_set(J,0,0,12*x1*x1 + 4*x2 - 42);
    gsl_matrix_set(J,0,1,4*x1 + 4*x2);
    gsl_matrix_set(J,1,0,4*x1 + 4*x2);
    gsl_matrix_set(J,1,1,12*x2*x2 + 4*x1 - 26);
	return GSL_SUCCESS;
}
