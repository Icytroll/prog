#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<math.h>
#include<time.h>
#include<string.h>

struct data {int N; int* sum; unsigned int seed;};

int randomDart(){
	double x = (double)rand()/(double)RAND_MAX-0.5;
	double y = (double)rand()/(double)RAND_MAX-0.5;
	double r = sqrt(x*x+y*y);
	return (r<=0.5 ? 1:0);
}

void* throwDarts(void* param){
	struct data data = *(struct data*)param;
	int N = data.N;
	int* sum = data.sum;
	unsigned int seed = data.seed;
	double x,y,r;
	for(int i=0;i<N;i++) {
		x = (double)rand_r(&seed)/(double)RAND_MAX-0.5;
		y = (double)rand_r(&seed)/(double)RAND_MAX-0.5;
		r = sqrt(x*x+y*y);
		if(r<=0.5) *sum+=1;
	}
	return NULL;
}

double Pthreads(int N, int n){
	
	int Nsum[n];
	for(int i=0;i<n;i++) Nsum[i]=0;
	int Ni = 0;
	struct data data[n];
	pthread_t thread[n];
	
	for(int i=0;i<n;i++){
		data[i].N = N/n;
		data[i].sum = &Nsum[i];
		data[i].seed = rand();
		pthread_create(&thread[i],NULL,throwDarts,(void*)&data[i]);
	}
	
	for(int i=0;i<n;i++){
		pthread_join(thread[i],NULL);
		Ni += Nsum[i];
	}
	
	return 4.0*Ni/N;
}

int main() {
	
	int n = 16; // threads
	int Nmax = 1e3, dN = 1e1, N;

	/* Pthreads */
	
	double Pi;
	printf("N      Pi\n");
	for(int i=0;i<Nmax/dN;i++){
		N = (i+1)*dN;
		Pi = Pthreads(N,n);
		printf("%.5i %.10g\n",N,Pi);
	}
	printf("\n\n");

	/* OpenMP */
	
	int Ni_openMP;
	for(int i=dN;i<=Nmax;i+=dN) {
		
		Ni_openMP = 0;
		#pragma omp parallel for reduction(+:Ni_openMP)
		for(int j=0;j<i;j++) {
			Ni_openMP += randomDart();
		}
		printf("%.5i %.10g\n",i,4.0*Ni_openMP/i);
	}
	
	return 0;
}
