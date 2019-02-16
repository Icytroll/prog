#include<stdio.h>
#include<stdlib.h>
#include"nvector.h"


nvector* nvector_alloc(int n){
	nvector* v = malloc(sizeof(nvector));
	(*v).size = n;
	(*v).data = malloc(n*sizeof(double));
	if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
	return v;
}

void nvector_free(nvector* v){ free(v->data); free(v);} /* v->data is identical to (*v).data */

void nvector_set(nvector* v, int i, double value){ (*v).data[i]=value; }

double nvector_get(nvector* v, int i){return (*v).data[i]; }

double nvector_dot_product(nvector* u, nvector* v){
	int n;
	for(int i=0;i<v->size;i++){
		n += v->data[i]+u->data[i];
	}
	return n;
}

int nvector_equal(nvector* a, nvector* b){
	int fail;
	double value;
	for(int i=0;i<a->size;i++){
		value = a->data[i] - b->data[i];
		if(value != 0){
			fail = 1;
			break;}
		else{fail = 0;}
		
	}
	if(fail){return 0;}
	else{return 1;}
}

void nvector_print(char* s, nvector* v) {
	printf("%s",s);
	printf("{");
	for(int i=0;i<v->size;i++){
		printf("%g ",v->data[i]);
	}
	printf("}\n");
}

void nvector_set_zero(nvector* v) {
	for(int i=0;i<v->size;i++){
		nvector_set(v,i,0);
	}
}

void nvector_add(nvector* a, nvector* b) {
	if(a->size!=b->size) {printf("a and b are not the same size!\n");}
	else{for(int i=0;i<a->size;i++){a->data[i] = a->data[i] + b->data[i];}}
}

void nvector_sub(nvector* a, nvector* b) {
	if(a->size!=b->size) {printf("a and b are not the same size!\n");}
	else{for(int i=0;i<a->size;i++){a->data[i] = a->data[i] - b->data[i];}}
}

void nvector_scale(nvector* a, double x) {
	for(int i=0;i<a->size;i++){a->data[i]*=x;}
}
