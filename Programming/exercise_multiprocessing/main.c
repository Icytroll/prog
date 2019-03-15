#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<math.h>
#include<time.h>
#include<string.h>

struct data {int N, sum; unsigned int seed;};

int randomDart(){
	double x = (double)rand()/(double)RAND_MAX-0.5;
	double y = (double)rand()/(double)RAND_MAX-0.5;
	double r = sqrt(x*x+y*y);
	return (r<=0.5 ? 1:0);
}

void* throwDarts(void* param){
	struct data* data = (struct data*)param;
	double x,y,r;
	for(int i=0;i<data->N;i++) {
		x = (double)rand_r(&(data->seed))/(double)RAND_MAX-0.5;
		y = (double)rand_r(&(data->seed))/(double)RAND_MAX-0.5;
		r = sqrt(x*x+y*y);
		if(r<=0.5) data->sum++;
	}
	return NULL;
}

double Pthreads(int N, int nthreads){
	
	int Ni = 0;
	struct data data[nthreads];
	pthread_t thread[nthreads];
	
	for(int i=0;i<nthreads;i++){
		data[i].N = N/nthreads;
		data[i].sum = 0;
		data[i].seed = rand();
		pthread_create(&thread[i],NULL,throwDarts,(void*)&data[i]);
	}
	
	for(int i=0;i<nthreads;i++){
		pthread_join(thread[i],NULL);
		Ni += data[i].sum;
	}
	
	return 4.0*Ni/N;
}

int main(int argc, char** argv) {
	int a=argc>1?atoi(argv[1]):11;
	int b=argc>2?atoi(argv[2]):19;
	
	int nthreads = 16; // threads
	int N;

	/* Pthreads */
	
	double Pi;
	printf("N      Pi\n");
	for(int i=a;i<b;i++){
		N = pow(2,i);
		Pi = Pthreads(N,nthreads);
		printf("%i %g\n",N,Pi);
	}
	printf("\n\n");

	/* OpenMP */
	
	int Ni_openMP=0;
	for(int iscale=a;iscale<b;iscale++){
		N = pow(2,iscale);
	
		struct data data[nthreads];
		for(int i=0;i<nthreads;i++){
			data[i].N = N/nthreads;
			data[i].sum = 0;
			data[i].seed = rand();
		}

		#pragma omp parallel for
		for(int i=0;i<nthreads;i++)
		{
			throwDarts((void*)&(data[i]));
		}
		
		int Ntot=0;
		for(int i=0;i<nthreads;i++){
			Ni_openMP += data[i].sum;
			Ntot += data[i].N;
		}
	
		printf("%i %g\n",N,4.0*Ni_openMP/Ntot/2);
	}
	
	return 0;
}


