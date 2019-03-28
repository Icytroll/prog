#include<math.h>
#include<gsl/gsl_matrix.h>

void matrix_print(gsl_matrix* A, char* s, FILE* stream);

int jacobi_single_row(gsl_matrix* D, int p) {
	int n_rot = 0, iter = 0, maxiter = 100, n = D->size1, flag = 1;
	double phi,c,s,D_pp,D_pq,D_qq,D_pi,D_qi,dpp,dqq;
	while (flag && iter<maxiter) {
		iter++;
		flag = 0;
		for(int q=p+1;q<n;q++) {
			D_pp = gsl_matrix_get(D,p,p);
			D_pq = gsl_matrix_get(D,p,q);
			D_qq = gsl_matrix_get(D,q,q);
			phi = atan2(2*D_pq,D_qq-D_pp)/2;
			c = cos(phi);
			s = sin(phi);
			for(int i=p;i<n;i++) {
				if(i!=p && i!=q) {
					D_pi = gsl_matrix_get(D,p,i);
					D_qi = gsl_matrix_get(D,q,i);
					gsl_matrix_set(D,p,i,c*D_pi-s*D_qi);
					gsl_matrix_set(D,i,p,c*D_pi-s*D_qi);
					gsl_matrix_set(D,q,i,s*D_pi+c*D_qi);
					gsl_matrix_set(D,i,q,s*D_pi+c*D_qi);
				}
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
	}
	return n_rot;
}
