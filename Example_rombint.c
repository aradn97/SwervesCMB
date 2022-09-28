#include<stdio.h>
#include<math.h>

#include "numericx.h"

double f(double x) {
    return cos(x);
}

int main() {
    numericx Numericx;
    double tol = 1e-8;
    double ans = Numericx.rombint(f, 0, 3.14159265359, tol);
    printf("%2.9e",ans);
}