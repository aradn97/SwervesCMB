#ifndef	_NEUTRINO_H_
#define	_NEUTRINO_H_


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "initval.h"
#include "variables.h"
#include "numericx.h"

double rl[nrhopn+1],pl[nrhopn+1];				// Logarithm of massive neutrino density, pressure 				
double drl[nrhopn+1],dpl[nrhopn+1];			   // Derivative of density and pressure					
double ddrl[nrhopn+1];								// Double derivatve 
double dlna,amin,amnu;

/***********************************************************************************************************************/
/*	 Compute massive neutrino density and pressure in units of the mean density of one flavor of massless neutrinos.	  */
/*  Use cubic splines to interpolate from a table.																							  */
/***********************************************************************************************************************/

void nu1(double a,double rrnu[])
{
	double rhonu,pnu,d;
//	int nrhopn=10000;																						// Tempurary variable			
      
   if(amnu==0.0)																							// No massive neutrinos present
   {
	   rhonu=1.0;
      pnu=1.0/3.0;
      }
      
	d=log(a/amin)/dlna+1.0;																						
   int i=(int)d;																							// Calculate floor of d	
   d=d-i;																									// Calculate fractional part of d
      
   if(i<0) 
   {																											// Use linear interpolation, bounded by results for massless neutrinos.
      rhonu=rl[1]+(d-1)*drl[1];
      pnu=pl[1]+(d-1)*dpl[1];
      rhonu=min(exp(rhonu),1.0);
      pnu=min(exp(pnu),0.3333333333);
		}
   else if(i>nrhopn)
   {	   																									// This should not happen, unless the user evolves to z<0!
      rhonu=rl[nrhopn-1]+(d+i-nrhopn)*drl[nrhopn-1];
      pnu=pl[nrhopn-1]+(d+i-nrhopn)*dpl[nrhopn-1];
      rhonu=exp(rhonu);
      pnu=exp(pnu);
      }
   else
   {	   																									// Cubic spline interpolation.
      rhonu=rl[i-1]+d*(drl[i-1]+d*(3.0*(rl[i]-rl[i-1])-2.0*drl[i-1]-drl[i]+d*(drl[i-1]+drl[i]+2.0*(rl[i-1]-rl[i]))));
      pnu=pl[i-1]+d*(dpl[i-1]+d*(3.0*(pl[i]-pl[i-1])-2.0*dpl[i-1]-dpl[i]+d*(dpl[i-1]+dpl[i]+2.0*(pl[i-1]-pl[i]))));
      rhonu=exp(rhonu);
      pnu=exp(pnu);
      }

	rrnu[0]=rhonu;
	rrnu[1]=pnu;
	}

	
/************************************************************************************************************************************/  
/*  Compute the perturbations of density, energy flux, pressure, and shear stress of one flavor of massive neutrinos, 					*/
/*  in units of the mean density of one flavor of massless neutrinos, by integrating over momentum.											*/
/************************************************************************************************************************************/  

void nu2(double a,double xsanpass[],double psi0[],double psi1[],double psi2[])
{
	double drhonu,fnu,dpnu,shearnu;
	double g0[4],g1[nqmax0+2],g2[nqmax0+2],g3[nqmax0+2],g4[nqmax0+2];
	
	double qmax=nqmax0-0.5;
  	double const2=5.68219698;																		// 7*(pi^4)/120.
	double gf1,gf2,gf3,gf4;
	double q,aq,v;
	
	if (nqmax==0)
	{
		drhonu=0.0;
      fnu=0.0;
      dpnu=0.0;
      shearnu=0.0;
      }
      
   else																									//  q is the comoving momentum in units of k_B*T_nu0/c.
   {
   	g1[1]=0.0;
      g2[1]=0.0;
      g3[1]=0.0;
      g4[1]=0.0;

      for(int iq=2;iq<=(nqmax0+1);iq++)
      {
   	   q=iq-1.5;
         aq=a*amnu/q;
         v=1.0/sqrt(1.0+aq*aq);
         g1[iq]=qdn[iq-2]*psi0[iq-2]/v;
         g2[iq]=qdn[iq-2]*psi0[iq-2]*v;
         g3[iq]=qdn[iq-2]*psi1[iq-2];
         g4[iq]=qdn[iq-2]*psi2[iq-2]*v;
   		}
   	
   	g0[0] = splint(g1,nqmax0+1);
      g0[1] = splint(g2,nqmax0+1);
      g0[2] = splint(g3,nqmax0+1);
      g0[3] = splint(g4,nqmax0+1);
      
      gf1=g1[nqmax0+1];
      gf2=g2[nqmax0+1];
      gf3=g3[nqmax0+1];
      gf4=g4[nqmax0+1];
      
      drhonu=(g0[0]+gf1*2.0/qmax)/const2;
      fnu=(g0[2]+gf3*2.0/qmax)/const2;
      dpnu=(g0[1]+gf2*2.0/qmax)/(const2*3.0);
      shearnu=(2.0/3.0)*(g0[3]+gf4*2.0/qmax)/const2;
   	}
   
	xsanpass[0] = drhonu;
	xsanpass[1] = fnu;
	xsanpass[2] = dpnu;
	xsanpass[3] = shearnu;
	}
	
/*****************************************************************************************************/
/*  Compute the time derivative of the mean density in massive neutrinos and the shear perturbation. */
/*****************************************************************************************************/

void nuder(double a,double rhonu,double xsanpass[],double adotoa,double psi2[],double psi2dot[])
{
	double rhonudot,shearnudot;
	double g1[nqmax0+2],dtemp;
	double q,aq,aqdot,v,vdot,g0,gf1;
	double qmax=nqmax0-0.5;
	int i;
	
	if (nqmax==0)
	{
		rhonudot=0.0;
      shearnudot=0.0;
      xsanpass[0]=rhonudot;
		xsanpass[2]=shearnudot;
      return;
     	}
		
	g1[1]=0.0;																							// q is the comoving momentum in units of k_B*T_nu0/c.
   for(int iq=2;iq<=(nqmax0+1);iq++)
   {
	   q=iq-1.5;
      aq=a*amnu/q;
      aqdot=aq*adotoa;
      v=1.0/sqrt(1.0+aq*aq);
      vdot=-aq*aqdot/pow((1.0+aq*aq),1.5);
      g1[iq]=qdn[iq-1]*(psi2dot[iq-1]*v+psi2[iq-1]*vdot);
      }
		
	g0 = splint(g1,nqmax0+1);
   gf1=g1[nqmax0+1];
   shearnudot=((g0+gf1*2.0/qmax)/constan)*(2.0/3.0);
	dtemp=log(a/amin)/dlna+1.0;
   i=int(dtemp);
   dtemp=dtemp-i;
   
   if (i<1)
   {
	   rhonudot=drl[1]+(dtemp-1)*ddrl[1];  															// Use linear interpolation
      }
   else if (i>nrhopn) 
   {
	   rhonudot=drl[nrhopn]+(dtemp+i-nrhopn)*ddrl[nrhopn];										// This should not happen, unless the user evolves to z<0!	
      }
   else
   {																										// Cubic spline interpolation for rhonudot.
      rhonudot=drl[i]+dtemp*(ddrl[i]+dtemp*(3.0*(drl[i+1]-drl[i])-2.0*ddrl[i]-ddrl[i+1]+dtemp*(ddrl[i]+ddrl[i+1]+2.0*(drl[i]-drl[i+1]))));
      }

	rhonudot=rhonu*adotoa*rhonudot/dlna;
	xsanpass[0]=rhonudot;
	xsanpass[2]=shearnudot;
   }
	

#endif	