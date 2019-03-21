#ifndef HAVE_QSPLINE_H /* this is necessary for multiple includes */

struct qspline {int n; double *x; double *y; double *a; double *b; double *c;};
typedef struct qspline qspline;

qspline* qspline_alloc(int n, double *x, double *y); /* allocates and builds the quadratic spline */
double qspline_eval(struct qspline *s, double z);        /* evaluates the prebuilt spline at point z */
double qspline_deriv(struct qspline *s, double z); /* evaluates the derivative of the prebuilt spline at point z */
double qspline_integ(struct qspline *s, double z);  /* evaluates the integral of the prebuilt spline from x[0] to z */
void qspline_free(struct qspline *s); /* free memory allocated in qspline_alloc */

#define HAVE_QSPLINE_H
#endif
