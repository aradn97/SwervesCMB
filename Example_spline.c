#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;
    double x[1001];
    double y[1001];
    double y2[1001];

    FILE *fp;
    /* creating a file to write the results to */
    fp =fopen("test.d","w");
    
    for (int i = 1; i <= 1000; i++) 
    {
        /* populating values of sin(x) from 0 to 2*pi into an array of 1000 values */
        /* h = 2*3.14/1000 is the step size */
        x[i] = 2*3.141592653*i/1000;
        y[i] = sin(x[i]);
    }

    Numericx.spline(x, y, 1000, 1.0, 1.0, y2);
        
    for (int i = 1; i <= 1000; i++) 
    {
        /* printing value of d^2(y[i])/d(y[i])^2 to the file test.d */
        fprintf(fp,"%lf %lf\n", y[i], y2[i]);
    }
    fclose(fp);
}