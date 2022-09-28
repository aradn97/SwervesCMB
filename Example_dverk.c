#include<stdio.h>
#include<math.h>

#include "numericx.h"

// function containing the differential equations
void func(int n,double x,double y[],double yprime[])
{
    yprime[1] = -0.4*y[1] - 16*y[2];
    yprime[2] = y[1];
}

int main() {
    numericx Numericx;
    double y[3] = {0.0, 0.5, 0.5};  // Initial conditions, note that y[0] is not used
    
    double c[25]={0.0};
    double **w; 
    double start, end;
    int indicator;



    indicator = 1;
   	FILE *fp;
	fp = fopen("dverk.d", "w");
    for(int i=0;i<1000;i++)
    {
        start=i*0.01;
        end = (i+1)*0.01;
        Numericx.dverk(3, func, start, y, end, 1e-9, indicator, c, 3, w);
        indicator = 3;
        fprintf(fp,"%f %f %e %e\n",start,end,y[1],y[2]);
    }
    fclose(fp);
}


