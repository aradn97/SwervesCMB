#include<stdio.h>
#include<math.h>
#include "doubleinflation.h"
#include "initval.h"

int main()
{
	FILE *fpcreatedata;
	double klotwo,khitwo,ksteptwo;
	double ktwo;

	klotwo = 1.0e-5;
	khitwo = 1.0;
	ksteptwo = log(khitwo/klotwo)/10000;
	fpcreatedata = fopen("powers.d","w");
	

	for(int i=0;i<10000;i++)
	{
		ktwo =  klotwo*exp(ksteptwo*i);
		//		interpoletepowers(ktwo);
		fprintf(fpcreatedata,"%e %e %e %e\n",ktwo,poweradiabatic(ktwo),powerisotherm(ktwo),powercross(ktwo));		
		}

}	