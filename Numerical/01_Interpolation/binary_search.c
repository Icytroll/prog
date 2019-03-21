#include<stdio.h>

int binary_search(int n, double *x, double z) {
	printf("z=%g\n",z);
	for(int i=0;i<n;i++)
		printf("x[%i]=%g\n",i,x[i]);
	int i=0,j=n-1,mid;
	while(j-i>1){
		mid = (j+i)/2;
		if(z>x[mid]) i = mid;
		else j = mid;
	}
	return i;
}
