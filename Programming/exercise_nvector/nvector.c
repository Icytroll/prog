#include<stdio.h>
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

double nvector_dot_product(nvector* u, nvector* n); {
	int n;
	for(i=0;i<v->size;i++){
		n += v->data[i]+u->data[i];
	}
	return n;
}

void nvector_equal(nvector* a, nvector* b){
	int fail;
	for(i=0;i<a->size;i++){
		if(a->data[i]-b->data[i] != 0){
			fail = 1;
			break;
		}
		
	}
	if(fail){return 0;}
	else{return 1;}
}
