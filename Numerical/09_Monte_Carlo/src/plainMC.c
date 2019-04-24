#include<math.h>
#include<stdlib.h>
#include<vector.h>
#include<matrix.h>
#define RND ((double)rand()/RAND_MAX)

void randomx(vector* x, vector* a, vector* b) {
	int n = x->size;
	double ai, bi;
	for(int i=0;i<n;i++) {
		ai = vector_get(a,i);
		bi = vector_get(b,i);
		vector_set(x,i,ai+RND*(bi-ai));
	}
}

void plainMC(
	double f(vector* x),
	vector* a,
	vector* b,
	int N,
	double* result,
	double* error) {
	
	int n = a->size;
	vector* x = vector_alloc(n);
	
	double V = 1;
	for(int i=0;i<n;i++) V *= vector_get(b,i)-vector_get(a,i);
	printf("V = %g\n",V);
	
	double sum = 0, sum2 = 0, fx;
	for(int i=0;i<N;i++) {
		randomx(x,a,b);
		fx = f(x);
		sum  += fx;
		sum2 += fx*fx;
	}
	double avr = sum/N, var = sum2/N-avr*avr;
	*result = avr*V;
	*error = sqrt(var/N)*V;

	free(x);
}
