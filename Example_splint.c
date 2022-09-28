#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;
    double y[1001];
    
    for (int i = 1; i <= 1000; i++) 
    {
        /* populating values of sin(x) from 0 to pi into an array of 1000 values */
        /* h = 3.14/1000 is the step size */
        y[i] = sin(3.141592653*i/1000);
    }
    printf("%lf", Numericx.splint(y, 1000)*(3.141592653/1000));
}