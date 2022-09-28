#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;
    double y[1001];
    double dy[1001];
    FILE *fp;
    /* creating a file to write the results to */
    fp =fopen("test.d","w");
    
    for (int i = 1; i <= 1000; i++) 
    {
        /* populating values of sin(x) from 0 to 2*pi into an array of 1000 values */
        /* h = 2*3.14/1000 is the step size */
        y[i] = sin(2*3.14*i/1000);
    }
    Numericx.splini(1000);
    Numericx.splder(y, dy, 1000);
        
    for (int i = 1; i <= 1000; i++) 
    {
        /* printing value of dy[i]/h to the file test.d */
        fprintf(fp,"%lf %lf\n", y[i], dy[i]/(2*3.14/1000));
    }
    fclose(fp);
}