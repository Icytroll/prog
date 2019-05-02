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
	int N = 100;
	vector* xi = vector_alloc(N);
	vector* yi = vector_alloc(N);
	for(int i=0;i<N;i++) {
		vector_set(xi,i,i/10.);
		vector_set(yi,i,sin(i/10.));
	}
	
	// A - Train network to minimze errors for 1D data
	
	int n = 3;
	double f(double x) {return x*exp(-x*x);}
	ann* network = ann_alloc(n,f);
	vector_set_all(network->data,1);
	/*vector_set(network->data,0,0);
	vector_set(network->data,2,0.5);
	vector_set(network->data,3,0);
	vector_set(network->data,5,0.5);
	vector_set(network->data,6,0);
	vector_set(network->data,8,0.5);
	*/printf("n = %d, f(1) = %g\n",network->n,network->f(1));
	ann_feed_forward(network,1);
	
	ann_train(network,xi,yi);

	ann_free(network);

	// B - Include derivative and anti-derivative
	
	
		
	// C - ??
	
	
	return 0;
}
