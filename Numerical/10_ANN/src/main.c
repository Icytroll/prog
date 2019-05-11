#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<unistd.h>
#include<assert.h>
#include"ann.h"
#include"vector.h"
#include"matrix.h"

/*--- FUNCTION DECLARATIONS ---*/

void broyden(
	double f(vector* x),
	void set_df(vector* x, vector* df),
	vector* x,
	double epsilon);

double CC_init(double f(double x), double a, double b, double acc, double eps);

/*--- MAIN PROGRAM ---*/

int main() {
	
	// generate tabulated data of the error function
	double xmin=-3;
	double xmax=3;
	int N = 201;
	vector* xi = vector_alloc(N);
	vector* yi = vector_alloc(N);
	for(int i=0;i<N;i++) {
		double x=xmin+i*(xmax-xmin)/(N-1);
		double y=erf(x);
		vector_set(xi,i,x);
		vector_set(yi,i,y);
	}

	// A - Train network to minimze errors for 1D data
	// B - Include derivative and anti-derivative

	int n = 30;
	double f(double x) {return x*exp(-x*x);}
	double df(double x) {return (1-2*x*x)*exp(-x*x);}
	
	// allocate network
	ann* network = ann_alloc(n,f,df);
	
	// distribute ai parameters across the x-axis
	vector_set_all(network->data,1);
	for(int i=0;i<n;i++){
		double a=xmin+(i+0.5)*(xmax-xmin)/n;
		vector_set(network->data,i*3,a);
	}
	
	// train the network
	printf("Optimizing network with numerical broyden's update:\n");
	vector_print(network->data,"p_start =",stdout);
	ann_train(network,xi,yi);
	vector_print(network->data,"p_final =",stdout);
	
	// declare function for Clenshaw-Curtis adaptive integration
	double integ_f(double x) {
		return ann_feed_forward(network,x);
	}

	// print values for plotting
	double x;
	for(int i=0;i<N;i++) {
		x = vector_get(xi,i);
		fprintf(stderr,"%g %g %g %g %g %g %g\n",
			x,
			vector_get(yi,i),ann_feed_forward(network,x),
			(2*exp(-x*x))/(sqrt(M_PI)),ann_deriv(network,x),
			x*erf(x)+exp(-x*x)/(sqrt(M_PI)),CC_init(integ_f,xmin,x,1e-6,1e-6));
	}
	
	// clean up
	ann_free(network);
	vector_free(xi);
	vector_free(yi);
	
	return 0;
}
