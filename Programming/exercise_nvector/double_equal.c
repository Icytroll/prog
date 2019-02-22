#include<math.h>
#include<stdlib.h>
#include<stdio.h>

int double_equal(double a, double b, double tau, double epsilon) {

    int res;
    if (abs(a-b) < tau) {res = 1;}
    else if (abs(a-b)/(abs(a)+abs(b)) < epsilon/2) {res = 1;}
    else {res = 0;}

    return res;
}
