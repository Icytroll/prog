#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrix.h"
#include"vector.h"

double F(vector* c, double funs(int,double), double x) {
	int m = c->size;
	double F = 0;
	for(int i=0;i<m;i++) F += vector_get(c,i)*funs(i,x);
	return F;
}

double funs(int i, double x) {
	switch(i) {
		case 0: return log(x);    break;
		case 1: return 1.0;    break;
		case 2: return x;      break;
		default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
	}
}

void lsfit(
	int n, double funs(int,double),
	double* x, double* y, double* dy,
	vector* c, matrix* S, vector* dc);

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

	// Calculate least square coefficients and their uncertainties via the covariance matrix
	int n = n_lines;
	vector* c = vector_alloc(m);
	matrix* S = matrix_alloc(m,m);
	vector* dc = vector_alloc(m);
	lsfit(n,funs,x,y,dy,c,S,dc);
	vector_print(c,"c =",stdout);
	matrix_print(S,"S =",stdout);
	vector_print(dc,"dc =",stdout);
	
	free(x);
	free(y);
	free(dy);
	
	// Write least squares function data to file 'A.txt'
	vector* c_lb = vector_alloc(m);
	vector* c_ub = vector_alloc(m);
	vector_sub(c,dc,c_lb);
	vector_add(c,dc,c_ub);
	
	double x1 = 0.1, xn = 9.9, dx = 0.1;
	
	FILE* Astream = fopen("A.txt","w");
	for(double x=x1;x<=xn;x+=dx) {
		fprintf(Astream,"%g %g %g %g\n",x,F(c,funs,x),F(c_lb,funs,x),F(c_ub,funs,x));
	}
	//fclose(Astream);


/*--- B Uncertainties of the fitting coefficients ---*/
	
	/*
	C = R'*R
	C = W*T
	Tx = W'*ei
	*/

	/*
	matrix* RT = matrix_alloc(m,m);
	matrix_transpose(R,RT);
	matrix_print(RT,"R' =",stdout);
	matrix* RTR = matrix_alloc(m,m);
	matrix_mult(RT,R,RTR);
	matrix_print(RTR,"R'*R =",stdout);
	
	matrix* invRTR = matrix_alloc(m,m);
	inverse(RTR,invRTR);
	matrix_print(invRTR,"inv(R'*R) =",stdout);
	*/
	

	/*
	matrix* RR = matrix_alloc(m,m);
	qr_gs_decomp(RTR,RR);
	matrix_print(RTR,"Q =",stdout);
	matrix_print(RR,"R =",stdout);
	qr_gs_inverse(RTR,RR,invRTR);
	matrix_print(invRTR,"inv(R'*R) =",stdout);
	matrix* invRTRRTR = matrix_alloc(m,m);
	matrix_mult(invRTR,RTR,invRTRRTR);	
	matrix_print(invRTRRTR,"inv(R'*R)*(R'*R) =",stdout);	
	*/
	
/*--- C Ordinary least-squares solution by thin singular-value decomposition ---*/
		

	
	//matrix_free(A);
	//matrix_free(Q);
	//matrix_free(R);
	//vector_free(b);
	vector_free(c);
	matrix_free(S);
	vector_free(dc);
	vector_free(c_lb);
	vector_free(c_ub);
	//vector_free(Ac);
	//matrix_free(QR);
	//matrix_free(RT);
	//matrix_free(RTR);
	//matrix_free(invRTR);
	return 0;
}
