#include<stdlib.h>
#include<math.h>
#include<float.h>
#include"ann.h"
#include"vector.h"


void broyden_numerical(double f(vector* x), vector* x, double dx, double epsilon);

ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double)) {
	ann* network = (ann*)malloc(sizeof(ann));
	network->n = number_of_hidden_neurons;
	network->f = activation_function;
	vector* data = vector_alloc(3*number_of_hidden_neurons);
	network->data = data;
	return network;
}

double ann_feed_forward(ann* network, double x) {
	int n = network->n;
	double sum = 0, ai, bi ,wi;
	for(int i=0;i<n;i++) {
		ai = vector_get(network->data,i*n);
		bi = vector_get(network->data,i*n+1);
		wi = vector_get(network->data,i*n+2);
		sum += network->f((x-ai)/bi)*wi;
	}
	return sum;
}

void ann_train(ann* network, vector* xi, vector* yi) {
	int N = xi->size;
	
	double ann_error(vector* x) {
		double sum = 0;
		for(int i=0;i<N;i++) {
			sum += pow(ann_feed_forward(network,vector_get(xi,i)) - vector_get(yi,i),2);
		}
		return sum/N;
	}
	
	double epsilon = 1e-6, dx = sqrt(DBL_EPSILON);
	broyden_numerical(ann_error,network->data,dx,epsilon);
}

void ann_free(ann* network) {
	vector_free(network->data);
	free(network);
}
