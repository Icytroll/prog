#ifndef HAVE_CSPLINE_H /* this is necessary for multiple includes */

struct cspline {int n; double *x; double *y; double *a; double *b; double *c; double *d;};
typedef struct cspline cspline;

cspline* cspline_alloc(int n, double *x, double *y); /* allocates and builds the quadratic spline */
double cspline_eval(struct cspline *s, double z);        /* evaluates the prebuilt spline at point z */
double cspline_deriv(struct cspline *s, double z); /* evaluates the derivative of the prebuilt spline at point z */
double cspline_integ(struct cspline *s, double z);  /* evaluates the integral of the prebuilt spline from x[0] to z */
void cspline_free(struct cspline *s); /* free memory allocated in qspline_alloc */

#define HAVE_CSPLINE_H
#endif
