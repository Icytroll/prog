#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

void matrix_print(gsl_matrix* A, char* s, FILE* stream);

int jacobi_classic(gsl_matrix* V, gsl_matrix* D) {
	int q, n_rot = 0, iter = 0, maxiter = 100, n = D->size1, flag = 1;
	double phi,c,s,D_pp,D_pq,D_qq,D_pi,D_qi,V_ip,V_iq,dpp,dqq;
	
	gsl_vector* large_index = gsl_vector_alloc(n-1);
	gsl_vector_set_all(large_index,n-1);
	for(int i=0;i<(n-1);i++) {
		for(int j=i+1;j<n;j++) {
			if(gsl_matrix_get(D,i,j)>gsl_matrix_get(D,i,gsl_vector_get(large_index,i))) {
				gsl_vector_set(large_index,i,j);
			}
		}
	}

	while (flag && iter<maxiter) {
		iter++;
		flag = 0;
		for(int p=0;p<(n-1);p++) {
			q = gsl_vector_get(large_index,p);
			D_pp = gsl_matrix_get(D,p,p);
			D_pq = gsl_matrix_get(D,p,q);
			D_qq = gsl_matrix_get(D,q,q);
			phi = atan2(2*D_pq,D_qq-D_pp)/2;
			c = cos(phi);
			s = sin(phi);
			for(int i=0;i<n;i++) {
				if(i!=p && i!=q) {
					D_pi = gsl_matrix_get(D,p,i);
					D_qi = gsl_matrix_get(D,q,i);
					gsl_matrix_set(D,p,i,c*D_pi-s*D_qi);
					gsl_matrix_set(D,i,p,c*D_pi-s*D_qi);
					gsl_matrix_set(D,q,i,s*D_pi+c*D_qi);
					gsl_matrix_set(D,i,q,s*D_pi+c*D_qi);
				}
				V_ip = gsl_matrix_get(V,i,p);
				V_iq = gsl_matrix_get(V,i,q);
				gsl_matrix_set(V,i,p,c*V_ip-s*V_iq);
				gsl_matrix_set(V,i,q,s*V_ip+c*V_iq);
			}
			dpp = c*c*D_pp-2*s*c*D_pq+s*s*D_qq;
			dqq = s*s*D_pp+2*s*c*D_pq+c*c*D_qq;
			if(dpp!=D_pp) {
				gsl_matrix_set(D,p,p,dpp);
				flag = 1;
			}
			if(dqq!=D_qq) {
				gsl_matrix_set(D,q,q,dqq);
				flag = 1;
			}
			gsl_matrix_set(D,p,q,0);
			gsl_matrix_set(D,q,p,0);
			n_rot++;
		}
		
		gsl_vector_set_all(large_index,n-1);
		for(int i=0;i<(n-1);i++) {
			for(int j=i+1;j<n;j++) {
				if(fabs(gsl_matrix_get(D,i,j))>fabs(gsl_matrix_get(D,i,gsl_vector_get(large_index,i)))) {
					gsl_vector_set(large_index,i,j);
				}
			}
		}
	}
	return n_rot;
}
