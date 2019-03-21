int binary_search(int n, double *x, double z);

void linear_param(int i,double *p,double *x,double *y) {
	p[1] = (y[i+1]-y[i])/(x[i+1]-x[i]);
	p[0] = y[i]-p[1]*x[i];
}

double linterp_integ(int n,double *x, double *y, double z) {
	
	int i = binary_search(n,x,z);
	double sum = 0;
	double p[2];

	int j = 0;
	while(j<=i) {
		linear_param(j,p,x,y);
		if(j==i) sum += (p[0]*z+0.5*p[1]*z*z)-(p[0]*x[j]+0.5*p[1]*x[j]*x[j]);
		else sum += (p[0]*x[j+1]+0.5*p[1]*x[j+1]*x[j+1])-(p[0]*x[j]+0.5*p[1]*x[j]*x[j]);
		j++;
	}

	return sum;
}
