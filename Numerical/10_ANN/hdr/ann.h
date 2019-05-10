#include"vector.h"

#ifndef HAVE_ANN_H /* this is necessary for multiple includes */
typedef struct {int n; double(*f)(double); double(*df)(double); vector* data;} ann;

ann* ann_alloc(int number_of_hidden_neurons, double(*f)(double), double(*df)(double));
double ann_feed_forward(ann* network, double x);
double ann_deriv(ann* network, double x);
void ann_train(ann* network, vector* xi, vector* yi);
void ann_free(ann* network);

#define HAVE_ANN_H
#endif
