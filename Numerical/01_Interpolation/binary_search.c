int binary_search(int n, double *x, double z) {
	int i=0,j=n-1,mid;
	while(j-i>1){
		mid = (j+i)/2;
		if(z>x[mid]) i = mid;
		else j = mid;
	}
	return i;
}
