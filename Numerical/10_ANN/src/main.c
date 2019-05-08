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
	double xmin=-2;
	double xmax=2;
	int N = 10;
	vector* xi = vector_alloc(N);
	vector* yi = vector_alloc(N);
	for(int i=0;i<N;i++) {
		double x=xmin+i*(xmax-xmin)/(N-1);;
		double y=exp(-x*x);
		vector_set(xi,i,x);
		vector_set(yi,i,y);
	}
	vector_print(xi,"xi =",stdout);
	vector_print(yi,"yi =",stdout);

	// A - Train network to minimze errors for 1D data
	
	int n = 1;
	double f(double x) {return exp(-x*x);}
	
	ann* network = ann_alloc(n,f);
	
	vector_set_all(network->data,1);
	
	for(int i=0;i<n;i++){
		double a=xmin+(i+0.5)*(xmax-xmin)/n+0.02;
		vector_set(network->data,i,a);
	}
	
	printf("n = %d, f(1) = %g\n",network->n,network->f(1));
	
	ann_feed_forward(network,1);
	
	ann_train(network,xi,yi);

	ann_free(network);

	// B - Include derivative and anti-derivative
	
	
		
	// C - ??
	
	
	return 0;
}
