#ifndef	_DARKENERGY_H_
#define	_DARKENERGY_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "variables.h"

int ndyn;													
double wdyn;											 

// ndyn
// 0 - if w = -1 and ndyn =1 or 3
// 1 - const eqn of state perturb
// 3 - const eqn of state unperturb
// 2 - w function of a perturb (tabulated)
// 4 - w function of a unperturb (tabulated)

// wdyn - w for dark energy


/******************************************************************************************************************/
/*		Calculate the dark energy density. 																									*/
/*		To be called by the function finitial(),cmbflat(),dtauda()																	*/	
/******************************************************************************************************************/
double dynrho(double a)
{
	double dynrho1;
	if(ndyn!=0)
	{
		if(ndyn==1 || ndyn==3)
		{
			dynrho1 = pow(a,(-3.0*(1.0 + wdyn)));
			}
		/*else																		//	Calculated from tabulated data. Skip it for this moment						
		{
			klo=(a-aw[1])/da+1;
         khi=klo+1;
         dlna=log(a/aw(klo));
         dynrho1=logrho[klo]-3.0*((1+wd[klo]-dwda[klo]*aw[klo])*dlna+dwda[klo]*(a-aw[klo]));
			dynrho1 = exp(dynrho1);
			}*/
		}		
	else
		dynrho1=1.0;
		
	return dynrho1;
	}


/************************************************************************************************************************************/
/*		Dark energy perturbation. The function will return w=p/rho for dark energy																		*/
/*		It will be used by finitial() 																																*/
/************************************************************************************************************************************/

double wdyn_func(double a)
{																					// equation of state for dark energy

	double wdyn_func1;
//	double klo,khi,as,bs;
	
	if(ndyn!=0)																	// usual Omega_Lambda ndyn=0
	{
		if(ndyn==1 || ndyn==3) 												// constant equation of state ndyn=1
			wdyn_func1 = wdyn;												// tables with a,w,dw/da  ndyn=2
		
	/*	else if(ndyn==2 || ndyn==4)										// Tables are provided -> ndyn=2,4
		{
			klo=(a-aw[1])/da+1;
         khi=klo+1;
         as=(aw[khi]-a)/da;
         bs=1.0-as;
         wdyn_func1=as*wd[klo]+bs*wd[khi];
			}*/
		}	
	else
		wdyn_func1=-1.0;
			
	return wdyn_func1;
	}


/************************************************************************************************************************************/
/*   Dark energy perturbation.																																		*/															
/*	  derivative of pressure,density,shear																															*/
/************************************************************************************************************************************/

double dyn_nrg(double a,double grho,double gpres,double delphi,double delpsi,double xsanpass[])
{
	double dgrho,dgprs,dgtheta;
	double grhophi,dphi,adot;
	
	grhophi = dynrho(a) * omegav * grhom;
   dphi = sqrt(grhophi)*a;
   adot = sqrt(grho/3.0);

   if(ndyn == 0)																							// No perturbation so all zero.
   {
	   dgrho = 0.0;
      dgprs = 0.0;
      dgtheta = 0.0;
      return 1;
      }
	
	if(ndyn==1 || ndyn==3)
	{
		dgrho = (1.0+wdyn)*dphi*(delpsi- 1.5*adot*(1.0-wdyn)*delphi);
      dgprs = (1.0+wdyn)*dphi*(delpsi+ 1.5*adot*(1.0-wdyn)*delphi);
      dgtheta =  (1.0+wdyn)*ak2*dphi*delphi;
      }
      
 /*else																										// If table is provided. Neglect for this moment
	{
		w=wdyn_func(a);
      dwda=dwda_func(a);
      if (w>-1.0)
      {
	      dgrho = (1.0+w)*dphi*(delpsi - 1.5*adot*(1.0-w)*delphi);
         dgprs = (1.0+w)*dphi*(delpsi + 1.5*adot*(1.0-w)*delphi) + dphi*dwda*adot*a*delphi;
         dgtheta =  (1.0+w)*ak2*dphi*delphi;
         }
		}*/
	
	xsanpass[0] = dgrho;
	xsanpass[1] = dgprs;
	xsanpass[2] = dgtheta;
	return 1;
	}

/************************************************************************************************************************************/
/* 	Calculate the derivatives on phi and psi due to the dark energy																					*/
/* 	This function will be used by the dunction fderivs()																									*/
/************************************************************************************************************************************/

void dyn_phi(double a,double hdot,double grho,double gpres,double phi,double psi,double ysanpass[])
{
	double grhophi,phid,sourceterm,adot,addot;
	double dphi,dpsi,coef1,coef2;
	
	grhophi = dynrho(a) * omegav * grhom;
   phid = sqrt(grhophi)*a;
   sourceterm = -0.5*hdot*phid;
   adot = sqrt(grho/3.0);

   addot = (grho - 3.0*gpres)/6.0;

   coef1 = 2.0*adot;

   if(ndyn==1 || ndyn==3)
 	   coef2 = -1.5*(1.0-wdyn)*(addot- adot*adot*(3.5+1.5*wdyn)) + ak2;
	/*else																						// Data from table This part is not done for this moment
	{
		wdyn=wdyn_func(a);
      dwda=dwda_func(a);
      
      if(wdyn > -1.0) 
      	coef1=coef1+dwda*adot*a/(wdyn+1.0);

      coef2 = -1.5*(1.0-wdyn)*(addot -adot*adot*(3.5+1.5*wdyn)) + ak2 + 3.0*(dwda*adot*a)*adot;
      }*/

	dphi = psi;
   
   if (wdyn > -1.0)
   	dpsi = -coef1*psi-coef2*phi+sourceterm;
   else
   	dpsi = 0.0;
   	
   ysanpass[0] = dphi;
   ysanpass[1] = dpsi;
   }
																				

#endif			