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

/*--- MAIN PROGRAM ---*/

int main() {
	
	// generate tabulated data of a sine function
	double xmin=-4;
	double xmax=4;
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
	
	int n = 30;
	double f(double x) {return x*exp(-x*x);}
	double df(double x) {return (1-2*x*x)*exp(-x*x);}
	
	ann* network = ann_alloc(n,f,df);
	vector_set_all(network->data,1);
	for(int i=0;i<n;i++){
		double a=xmin+(i+0.5)*(xmax-xmin)/n;
		vector_set(network->data,i*3,a);
	}
	
	printf("Optimizing network with numerical broyden's update:\n");
	vector_print(network->data,"p_start =",stdout);
	ann_train(network,xi,yi);
	vector_print(network->data,"p_final =",stdout);

	for(int i=0;i<N;i++) {
		fprintf(stderr,"%g %g %g %g\n",vector_get(xi,i),vector_get(yi,i),ann_feed_forward(network,vector_get(xi,i)),ann_deriv(network,vector_get(xi,i)));
	}
	
	ann_free(network);
	
	// B - Include derivative and anti-derivative
	
	
		
	// C - ??
	
	
	return 0;
}
