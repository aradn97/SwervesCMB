#include<stdio.h>
#include<math.h>

#include "lensing.h"

double f(double x) {
    return cos(x);
}

int main() {
    lensing Lensing;
    double tol = 1e-8;
    double ans = Lensing.rombint1(f, 0, 3.14159265359, tol);
    printf("%2.9e",ans);
}