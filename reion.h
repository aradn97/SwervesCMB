#ifndef	_REION_H_
#define	_REION_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "variables.h"
#include "others.h"

/**********************************************************************************************************************/
/*  This function calculates the redshift of reionization from an optdlss and sets the reionization fraction rif=1	 */
/*  This function will be used by cmbflat() for calculaion when reionization starts and when it stops						 */ 
/*	 Thif function will use the function dtauda()																							 */
/**********************************************************************************************************************/

void reiopar()
{
	double da,akthom,a,optd;
   zri=0.0;
   rif=0.0;
   
   if (optdlss!=0.0)
   {
	   akthom=(2.3048e-9)*(1-yhe)*omegab*h0*h0;
      rif=1.0;																// Calculating the redshift of reionation the parameter needed by CMBANUS.
      da=0.00001;
      optd=0.0;
      a=1-da;
		for(int na=1;optd<optdlss;na++)
		{
      	optd=optd+da*rif*akthom*dtauda(a)/(a*a);
			a=1-na*da;
			}
      zri=1.0/a-1.0;
      zristp=0.07*zri-1.0;												//stop time of reionization

     if (zristp<0.0) 
      	zristp=0.0;						 								//mininum redshift can be 0
      }
	}  

#endif
