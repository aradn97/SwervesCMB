#include <stdio.h>
#include <math.h>

#include "numericx.h"

double **matrix(int col, int row)
{
    double **m;
    m = (double **)malloc(col * sizeof(double));
    if (!m)
    {
        printf("allocation failure 1 in dmatrix()");
        exit(1);
    }

    /* allocate rows and set pointers to them */
    for (int i = 0; i <= col; i++)
    {
        m[i] = (double *)malloc(row * sizeof(double));
        if (!m[i])
        {
            printf("allocation failure 2 in dmatrix()");
            exit(1);
        }
    }
    printf("Successfully allocated\n");
    /* return pointer to array of pointers to rows */
    return m;
}

int main()
{
    numericx Numericx;
    double x[1001];
    double **y;
    double **y2;

    y = matrix(1000, 3);
    y2 = matrix(1000, 3);

    FILE *fp;
    /* creating a file to write the results to */
    fp = fopen("test_spline2D.d", "w");

    for (int i = 0; i < 1000; i++)
    {
        /* populating values of sin(x) from 0 to 2*pi into an array of 1000 values */
        /* h = 2*3.14/1000 is the step size */
        x[i] = 2 * 3.141592653 * i / 1000;
        y[i][0] = sin(x[i]);
        y[i][1] = sin(2 * x[i]);
        y[i][2] = sin(3 * x[i]);
    }

    Numericx.spline2D(x, y, 1000, 3, 1.0, 1.0, y2);

    for (int i = 0; i < 1000; i++)
    {
        /* printing value of d^2(y[i])/d(y[i])^2 to the file test.d */
        fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", x[i], y[i][0], y2[i][0], y[i][1], y2[i][1], y[i][2], y2[i][2]);
    }
    fclose(fp);
}