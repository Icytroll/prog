int binary_search(int n, double *x, double z);

double linterp(int n, double *x, double *y, double z) {
	
	int i = binary_search(n,x,z);
	double a = (y[i+1]-y[i])/(x[i+1]-x[i]);
	return y[i]+a*(z-x[i]);
	
}
