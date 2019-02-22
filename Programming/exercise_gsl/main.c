#include<stdio.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_airy.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

int main() {
	
	/* Make a data file for Ai and Bi */

	FILE * fid;
	double aa = -15, bb = 5, n = 1000;
	double dx = (bb-aa)/(n-1);

	fid = fopen("data.txt","w");
	for (int i = 0;i<n;i++) {
		fprintf(fid,"%g %g %g\n",aa+i*dx,gsl_sf_airy_Ai(aa+i*dx,GSL_PREC_DOUBLE),gsl_sf_airy_Bi(aa+i*dx,GSL_PREC_DOUBLE));
	}

	fclose(fid);


	/* Solve linear algebro problem */

	fprintf(stdout,"Linear algebra problem:\n\n");
	
	double A_data[] = { 6.13, -2.90,  5.86,
			            8.08, -6.31, -3.89,
			           -4.36,  1.00,  0.19};

	double b_data[] = { 6.23,  5.37,  2.29};

	gsl_matrix_view A = gsl_matrix_view_array (A_data,3,3);
	gsl_vector_view b = gsl_vector_view_array (b_data,3);
	gsl_vector *x = gsl_vector_alloc (3);

	int s;

	gsl_permutation * p = gsl_permutation_alloc (3);
	gsl_linalg_LU_decomp (&A.matrix, p, &s);
	gsl_linalg_LU_solve  (&A.matrix, p, &b.vector, x);
	
	fprintf(stdout,"x = \n");
	gsl_vector_fprintf(stdout, x, "%.3g");
	fprintf(stdout,"\n");
	
	// redefine A-matrix for normal matrix operations
	double A_sol[3][3] = {{ 6.13, -2.90,  5.86},
	                      { 8.08, -6.31, -3.89},
	                      {-4.36,  1.00,  0.19}};

	double b_sol[3] = {0};
	
	// manual multiplication
	fprintf(stdout,"Multiplication of A and x:\n");
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			b_sol[i] += A_sol[i][j]*gsl_vector_get(x,j);
		}
		fprintf(stdout,"%g\n",b_sol[i]);
	}
	
	fprintf(stdout,"Should be:\n");
	for(int i=0;i<3;i++)
		fprintf(stdout,"%g\n",b_data[i]);
	
	// free up memory
	gsl_permutation_free(p);
	gsl_vector_free(x);
	

	/* Eigenvalues and eigenvectors of 4th order Hilbert matrix */
	
	fprintf(stdout,"\n\nEigenvalues and eigenvectors of Hilbert matrix:\n\n");

	// allocate variables
	int N = 4;
	gsl_matrix* m    = gsl_matrix_alloc(N,N);
	gsl_matrix* evec = gsl_matrix_alloc(N,N);
	gsl_vector* eval = gsl_vector_alloc(N);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(N);

	// make hilbert matrix
	double val;
	fprintf(stdout,"Hilbert matrix:\n");
	for(double i = 0;i<N;i++) {
		for(double j = 0;j<N;j++) {
			val = 1/(i+j+1);
			gsl_matrix_set(m,i,j,val);
			fprintf(stdout,"%5.3f ",val);
		}
		fprintf(stdout,"\n");
	}

	// solve for eigenvalues and eigenvectors
	gsl_eigen_symmv(m,eval,evec,w);

	// sort the results
	gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
	
	// printf results to stdout
	fprintf(stdout,"\n");
	for(int i = 0;i<N;i++) {
		fprintf(stdout,"eigenvalue \n%g\n",gsl_vector_get(eval,i));
		fprintf(stdout,"eigenvector\n");
		gsl_vector_view evec_i = gsl_matrix_column(evec,i);
		gsl_vector_fprintf(stdout,&evec_i.vector,"%g");
		fprintf(stdout,"\n");
	}

	// free up memory	
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(m);
	gsl_matrix_free(evec);
	gsl_vector_free(eval);

	return 0;
}
