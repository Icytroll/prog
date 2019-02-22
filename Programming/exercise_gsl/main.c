#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_airy.h>
#include<gsl/gsl_linalg.h>

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

	printf("x = \n");
	gsl_vector_fprintf(stdout, x, "%.3g");
	printf("\n");
	
	// redefine A-matrix for normal matrix operations
	double A_sol[3][3] = {{ 6.13, -2.90,  5.86},
	                      { 8.08, -6.31, -3.89},
	                      {-4.36,  1.00,  0.19}};

	double b_sol[3] = {0};
	
	// Manual multiplication
	printf("Multiplication of A and x:\n");
	for(int i=0;i<3;i++) {
		for(int j=0;j<3;j++) {
			b_sol[i] += A_sol[i][j]*gsl_vector_get(x,j);
		}
		printf("%g\n",b_sol[i]);
	}
	
	printf("Should be:\n");
	for(int i=0;i<3;i++)
		printf("%g\n",b_data[i]);
	
	// free up memory
	gsl_permutation_free(p);
	gsl_vector_free(x);
	

	return 0;
}
