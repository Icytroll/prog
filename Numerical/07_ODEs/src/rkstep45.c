#include"vector.h"
#include"matrix.h"

void rkstep45(
	double t,
	double h,
	vector* yt,
	void f(double t, vector* y, vector* dydt, void* params),
	void* params,
	vector* yth,
	vector* err) {
	
	// Allocate temporary y-array and slope-arrays
	int n = yt->size;
	vector* y = vector_alloc(n);
	vector* k0 = vector_alloc(n);
	vector* k1 = vector_alloc(n);
	vector* k2 = vector_alloc(n);
	vector* k3 = vector_alloc(n);
	vector* k4 = vector_alloc(n);
	vector* k5 = vector_alloc(n);
	
	f(t,      yt,k0,params); 
	for(int i=0;i<n;i++) vector_set(y,i,vector_get(yt,i)+1./4*vector_get(k0,i)*h);
	f(t+1./2*h,y,k1,params); 
	for(int i=0;i<n;i++) {
		vector_set(y,i,vector_get(yt,i)+
		(3./32*vector_get(k0,i)
		+9./32*vector_get(k1,i))*h);
	}
	f(t+3./8*h,y,k2,params); 
	for(int i=0;i<n;i++) {
		vector_set(y,i,vector_get(yt,i)+
		(1932./2197*vector_get(k0,i)
		-7200./2197*vector_get(k1,i)
		+7296./2197*vector_get(k2,i))*h);
	}
	f(t+12./13*h,y,k3,params); 
	for(int i=0;i<n;i++) {
		vector_set(y,i,vector_get(yt,i)+
		(439./216*vector_get(k0,i)
		-8*vector_get(k1,i)
		+3680./513*vector_get(k2,i)
		-845./4104*vector_get(k3,i))*h);
	}
	f(t+h,y,k4,params); 
	for(int i=0;i<n;i++) {
		vector_set(y,i,vector_get(yt,i)+
		(-8./27*vector_get(k0,i)
		+2*vector_get(k1,i)
		-3544./2565*vector_get(k2,i)
		+1859./4014*vector_get(k3,i)
		-11./40*vector_get(k4,i))*h);
	}
	f(t+1./2*h,y,k5,params); 
	for(int i=0;i<n;i++) {
		vector_set(yth,i,vector_get(yt,i)+
		(16./135*vector_get(k0,i)
		+6656./12825*vector_get(k2,i)
		+28561./56430*vector_get(k3,i)
		-9./50*vector_get(k4,i)
		+2./55*vector_get(k5,i))*h);
		vector_set(y,i,vector_get(yt,i)+
		(25./216*vector_get(k0,i)
		+1408./2565*vector_get(k2,i)
		+2197./4104*vector_get(k3,i)
		-1./5*vector_get(k4,i))*h);
		vector_set(err,i,vector_get(yth,i)-vector_get(y,i));
	}

	vector_free(y);
	vector_free(k0);
	vector_free(k1);
	vector_free(k2);
	vector_free(k3);
	vector_free(k4);
	vector_free(k5);
}
