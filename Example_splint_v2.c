#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;
    double xa[1001];
    double ya[1001];
    double y2a[1001];
    
    for (int i = 1; i <= 1000; i++) 
    {
        xa[i] = 3.141592653*i/1000;
        /* populating values of sin(x) from 0 to pi into an array of 1000 values */
        /* h = 3.14/1000 is the step size */
        ya[i] = sin(xa[i]);
    }

    /* populating values the second derivative of sin(x) from 0 to pi into an array of 1000 values */
    Numericx.spline(xa, ya, 1000, 1.0, 1.0, y2a);

    /* Interpolate the value of sin(x) at x = pi*(0.5+0.001) */
    printf("%lf", Numericx.splint_v2(xa, ya, y2a, 1000, 3.141592653*501/1000));
}