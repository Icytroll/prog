#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrix.h"
#include"vector.h"

double dF(matrix* S, double funs(int,double), double x) {
	int m = S->size1;
	double dF = 0;
	for(int i=0;i<m;i++)
		for(int j=0;j<m;j++)
			dF += funs(i,x)*matrix_get(S,i,j)*funs(j,x);
	return sqrt(dF);
}

double F(vector* c, double funs(int,double), double x) {
	int m = c->size;
	double F = 0;
	for(int i=0;i<m;i++) F += vector_get(c,i)*funs(i,x);
	return F;
}

double funs1(int i, double x) {
	switch(i) {
		case 0: return log(x);    break;
		case 1: return 1.0;    break;
		case 2: return x;      break;
		default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
	}
}

double funs2(int i, double x) {
	switch(i) {
		case 0: return log(x);    break;
		case 1: return x;    break;
		case 2: return x*x;      break;
		default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
	}
}

void lsfit(
	int n, double funs(int,double),
	double* x, double* y, double* dy,
	vector* c, matrix* S);

void inverse(matrix* A, matrix* invA);



int main() {
	
/*--- A Ordinary least-squares fit by QR-decomposition ---*/

	// Number of fitting functions
	int m = 3;
	
	// Make temporary arrays (in case we somehow want to change the amount of data points,very unnecessary here, just trying it out)
	double* x = malloc(100*sizeof(double));
	double* y = malloc(100*sizeof(double));
	double* dy = malloc(100*sizeof(double));

	// Read data
	FILE* fid = fopen("input.data","r");
	int n_lines = 0;
	while (fscanf(fid,"%lf %lf %lf",&x[n_lines],&y[n_lines],&dy[n_lines]) == 3) {;
		n_lines++;
	}
	fclose(fid);

/*--- B Uncertainties of the fitting coefficients ---*/
	
	// Calculate least square coefficients and their uncertainties via the covariance matrix
	int n = n_lines;
	vector* c = vector_alloc(m);
	matrix* S = matrix_alloc(m,m);
	
	lsfit(n,funs1,x,y,dy,c,S);
	printf("First set of functions\n\n");
	vector_print(c,"Coefficients =",stdout);
	matrix_print(S,"Covariance matric =",stdout);
	
	// Write first set of functions data
	double x1 = 0.1, xn = 9.9, dx = 0.1;
	
	FILE* fitstream = fopen("fit.txt","w");
	for(double x=x1;x<=xn;x+=dx) {
		fprintf(fitstream,"%g %g %g %g\n",x,F(c,funs1,x),F(c,funs1,x)+dF(S,funs1,x),F(c,funs1,x)-dF(S,funs1,x));
	}
	
	
	lsfit(n,funs2,x,y,dy,c,S);
	printf("Second set of functions\n\n");
	vector_print(c,"Coefficients =",stdout);
	matrix_print(S,"Covariance matrix =",stdout);
	
	// Write second set of functions data
	fprintf(fitstream,"\n\n");
	for(double x=x1;x<=xn;x+=dx) {
		fprintf(fitstream,"%g %g %g %g\n",x,F(c,funs2,x),F(c,funs2,x)+dF(S,funs2,x),F(c,funs2,x)-dF(S,funs2,x));
	}

	fclose(fitstream);


	


/*--- C Ordinary least-squares solution by thin singular-value decomposition ---*/
		

	
	free(x);
	free(y);
	free(dy);
	vector_free(c);
	matrix_free(S);
	
	return 0;
}
