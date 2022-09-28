#include<stdio.h>
#include<math.h>

#include "numericx.h"


int main() {
    numericx Numericx;

    double x[10001];
    double y[10001];
    int n = 2;

    FILE *fp;
    /* creating a file to write the results to */
    fp = fopen("bessj.d","w");
    
    for (int i = 0; i <= 10000; i++) 
    {
        /* populating the x array 0 to 20 with h = 20/10000 step size */
        x[i] = 20.0*i/10000.0;
        /* calculate the Bessel function J_n(x) where 0 <= x <= 20 and n = 2 */
        y[i] = Numericx.bessj(n, x[i]);
    }
        
    for (int i = 0; i <= 10000; i++) 
    {
        /* printing the result of J_2(x) to the file bessj.d */
        fprintf(fp,"%e %e\n", x[i], y[i]);
    }
    fclose(fp);
}