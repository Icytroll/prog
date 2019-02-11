#include<stdio.h>
#include<stdlib.h>
#include"komplex.h"

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

komplex komplex_sub (komplex a, komplex b) {
	komplex results = {a.re - b.re , a.im - b.im};
	return results;
}

int komplex_equal (komplex a, komplex b, double acc, double eps) {
	if (abs(a.re-b.re) < acc && abs(a.im-b.im) < acc) {return 1;}
	else if (abs(a.re-b.re)/(abs(a.re)+abs(b.re)) < eps/2 && abs(a.im-b.im)/(abs(a.im)+abs(b.im)) < eps/2) {return 1;}
	else {return 0;}
}
