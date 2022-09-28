#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;
    double xa[1001];
    double ya[1001];
    int n = 1000;
    
    for (int i = 1; i <= 1000; i++) 
    {
        xa[i] = 3.141592653*i/1000;
        /* populating values of sin(x) from 0 to pi into an array of 1000 values */
        /* h = 3.14/1000 is the step size */
        ya[i] = sin(xa[i]);
    }

    /* muptiply by the step size h to find the integrand over the array (sin(x) from 0 to pi) */
    double z = (3.141592653/1000)*Numericx.integ(ya, n);

    /* print the value of the integrand */
    printf("%lf", z);
}