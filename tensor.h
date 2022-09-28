#ifndef TENSOR_H_
#define TENSOR_H_



#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "initval.h"
#include "variables.h"
#include "numericx.h"
#include "darkenergy.h"
#include "neutrino.h"
#include "others.h"
#include "reion.h"
#include "recombination.h"
//#include "recfast.h"
#include "lensing.h"


class tensor {

// double yhe = Variables.yhe;
// double h0 = Variables.h0;
// double omegab = Variables.omegab;
// double tcmb = Variables.tcmb;
// double omegavdyn = Variables.omegavdyn;
// double omegav = Variables.omegav;
// double omegac = Variables.omegac;
// double annur = Variables.annur;
// double annunr = Variables.annunr;
// double rif = Variables.rif;
// double zri = Variables.zri;
// double grhom = Variables.grhom;
// double grhonr = vpc.Variable_grhonr(a);
// double grhor = vpc.Variable_grhor(a);
// double grhog = vpc.Variable_grhog(a);
// double dimflag = Variables.dimflag;



// int ndyn = Variables.ndyn;
// int initfl = Variables.initfl;


// int lmaxg = Variables.lmaxg;
// int lmaxnr = Variables.lmaxnr;
// int lmaxnu = Variables.lmaxnu;
// int lmaxt = Variables.lmaxt;
// int iq0 = Variables.iq0;
// int iq1 = Variables.iq1;
// int iq2 = Variables.iq2;
// int nqmax = Variables.nqmax;

// double amnu = Variables.amnu;

public:    


/************************************************************************************************************************************/
/*		This function will initiallize the tensor perturbation terms																						*/
/************************************************************************************************************************************/

void finitialt(double tau,double y[])
{
	int ind1,ind2,ind3;
	double a;
//   printf("\nTHIS ");
   a=tau*Variables.adotrad;
   y[1]=a;

																										// Tensor modes
   y[2]=1.0;																							// ht=1.0
   y[3]=0.0;																							// htpr=0.0
  // if(fabs(Variables.omegak)<0.0001)
   {
      ind1=4;
      ind2=ind1+Variables.lmaxt+1;
      ind3=ind2+Variables.lmaxt+1;
      }
  // else   
  // {
  //    ind1=2;
  //    ind2=ind1+Variables.lmaxt-1;
  //    ind3=ind2+Variables.lmaxt-1;
   //   }

   for(int l=0;l<=Variables.lmaxt;l++)
   {
 //     printf("\n%d %d",ind1+l,ind1+2);
   	y[ind1+l]=0.0;
   	y[ind2+l]=0.0;
   	}
   for(int l=0;l<=lmaxnr0;l++)
   	y[ind3+l]=0.0;

//   printf("I AM HERE %e",y[2]);   exit(1);
}


/************************************************************************************************************************************/
/*  This function will calculate the time derivatives of the tensor perturbation                                                    */
/*  This part is based on the following papers :                                                                                   */
/*                                                                                                                                  */
/*	 1) Cosmology                                                                                                               */
/*   Steven Weinberg (Page No  310-330)                                                                                             */
/************************************************************************************************************************************/

void fderivst(int n,double x,double y[],double yprime[])
{
	double ep,tau,a,a2,tcp,tcp1,tcp2,rhonu,pnu,omegavdyn,grho;
   double adotoa;
	double shearg,shearr,pi,ht,htpr,psie,deltat0,deltap0;
	double ep0=1.0e-2,xsanpass[3],htdpr,opac1,ysanpass[2];
	int ind1,ind2,ind3;
   double sourcet,sourcee,sourceb,p2;
   double sqrt6 = 2.4494897;
   //double p2,sourcet,sourcee,sourceb;
   double ak, ak2;

   recombination Recombination;
   neutrino Neutrino;
   darkenergy Darkenergy;
   VPC vpc;

   curvature Curvature(-Variables.omegak);
   
   ak = Variables.ak;
   ak2 = Variables.ak2;
    
   if(fabs(Variables.omegak)>0.0001)
   {
      ak = sqrt(ak2 - 3*Variables.curv);
      ak2 = ak*ak;
   } 

   if (Variables.ak>0.06*Variables.epsw)																					// ep is used to stop the tight coupling approximation.
   	ep=ep0;
   else
   	ep=1.17*ep0;

	tau=x;
   a=y[1];																								// Time of the calculation

   a2=a*a;

   Recombination.thermo(tau,xsanpass);																			// Calculate Opacity
   opac1=xsanpass[1];
   double opac = opac1;
																											// Tight Coupling parameters
   tcp=0.0;
   tcp1=Variables.ak/opac1;
   tcp2=1.0/(opac1*tau);
   if ((tcp1>ep)||(tcp2>ep))
   	tcp=1.0;																							// Not tightly Coupled

																											// Compute expansion rate.
   if (Variables.amnu==0)																						// No massive Neutronos
   {
   	rhonu=1.0;
      pnu=1.0/3.0;
      }
   else
   	Neutrino.nu1(a,ysanpass);																				// Masssive Neutrinos Present, So calculate pressure and density
   rhonu = ysanpass[0];
   pnu = ysanpass[1];

   omegavdyn = Variables.omegav*Darkenergy.dynrho(a);

   if(Variables.VPCflag == 1)  
   {
      omegavdyn = omegavdyn * vpc.Variable_Lambda(a);
      }
   																											//  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
   grho=vpc.Variable_grhom(a)*(Variables.omegac+Variables.omegab)/a+(vpc.Variable_grhog(a)+vpc.Variable_grhor(a)*Variables.annur+vpc.Variable_grhonr(a)*Variables.annunr*rhonu)/a2+vpc.Variable_grhom(a)*omegavdyn*a2;

	if(Variables.dimflag==1)																						// Energy density for  5 dimensions
	{
		grho=pow((sqrt(grho)+sqrt(Variables.omegav*vpc.Variable_grhom(a)*a2)),2);
		}

   adotoa=sqrt(grho/3.0);
   yprime[1]=adotoa*a;


	/*********************************************************************************************************************************/
	/*  Relativistic components																																		*/
	/*********************************************************************************************************************************/


	/*********************************************************************************************************************************/
	/*  Tensors																																								*/
	/*********************************************************************************************************************************/

   ht=y[2];
   htpr=y[3];
   yprime[2]=htpr;
   if(fabs(Variables.omegak)<0.0001)
      htdpr=-2*adotoa*htpr-Variables.ak2*ht+24*pi;															//Eq(W.73), Eq(R.3)
   else
      htdpr=-2*adotoa*htpr - (ak2 + 2.0*Variables.curv)*ht;      

   yprime[3]=htdpr;



   if(fabs(Variables.omegak)<0.0001)
   {
      ind1=4;
      ind2=ind1+Variables.lmaxt+1;
      ind3=ind2+Variables.lmaxt+1;
      shearg=y[ind1]/15+y[ind1+2]/21+y[ind1+4]/35;
      shearr=y[ind3]/15+y[ind3+2]/21+y[ind3+4]/35;
      pi=(vpc.Variable_grhog(a)*shearg+Variables.annur*vpc.Variable_grhor(a)*shearr)/(a2*3);												//
      }
   else
   {
      ind1=2;
      ind2=ind1+Variables.lmaxt-1;
      ind3=ind2+Variables.lmaxt-1;
      p2=0.1*(y[ind1+2]-sqrt6*y[ind2+2]);
      sourcet=opac*p2-htpr;
      sourcee=-opac*sqrt6*p2;
      sourceb=0.0;
      }   

	/*********************************************************************************************************************************/
	/*  Photon perturbations																																			*/
	/*********************************************************************************************************************************/
	psie=y[ind1]/10.0+y[ind1+2]/7.0+3*y[ind1+4]/70.0-3.0*y[ind2]/5.0+6.0*y[ind2+2]/7.0-3.0*y[ind2+4]/70.0;
																											//Eq(Y.61), Eq(R.3)
   
   if(fabs(Variables.omegak)<0.0001){
   if (tcp==1) 																						// Not Tightly Couppled
   {
	   yprime[ind1]=-Variables.ak*y[ind1+1]-opac1*y[ind1]+opac1*psie-htpr;							// Eq(R.3)
      yprime[ind2]=-Variables.ak*y[ind2+1]-opac1*y[ind2]-opac1*psie;									// Eq(R.3)
      
      if(Variables.VPCflag == 1)
      {
         yprime[ind1] = yprime[ind1] - vpc.Variable_G_c(a) * y[ind1]; 
         yprime[ind2] = yprime[ind2] - vpc.Variable_G_c(a) * y[ind2];          
         }
      for(int l=1;l<=Variables.lmaxt-1;l++)
      {
	      yprime[ind1+l]=Variables.ak*(l*y[ind1-1+l]-(l+1)*y[ind1+1+l])/(2*l+1)-opac1*y[ind1+l]; // Ma, Bertschinger Eq(64)
         yprime[ind2+l]=Variables.ak*(l*y[ind2-1+l]-(l+1)*y[ind2+1+l])/(2*l+1)-opac1*y[ind2+l]; // Ma, Bertschinger Eq(64)

         if(Variables.VPCflag == 1)
         {
            yprime[ind1+l] = yprime[ind1+l] - vpc.Variable_G_c(a) * y[ind1+l]; 
            yprime[ind2+l] = yprime[ind2+l] - vpc.Variable_G_c(a) * y[ind2+l];          
            }
         }

																											// Truncate moment expansion
																											// Ma, Bertschinger Eq(64)
      yprime[ind1+Variables.lmaxt]=Variables.ak*y[ind1-1+Variables.lmaxt]-(Variables.lmaxt+1)/tau*y[ind1+Variables.lmaxt]-opac1*y[ind1+Variables.lmaxt];
      yprime[ind2+Variables.lmaxt]=Variables.ak*y[ind2-1+Variables.lmaxt]-(Variables.lmaxt+1)/tau*y[ind2+Variables.lmaxt]-opac1*y[ind2+Variables.lmaxt];
      if(Variables.VPCflag == 1)
      {
         yprime[ind1+Variables.lmaxt] = yprime[ind1+Variables.lmaxt] - vpc.Variable_G_c(a) * y[ind1+Variables.lmaxt]; 
         yprime[ind2+Variables.lmaxt] = yprime[ind2+Variables.lmaxt] - vpc.Variable_G_c(a) * y[ind2+Variables.lmaxt];          
         }
      }

	else																									// Tightly Couppled
   {
	   deltat0=-4.0*htpr/(opac1*3.0);
      deltap0=-deltat0/4.0;
      y[ind1]=deltat0;
      y[ind2]=deltap0;

      if(Variables.VPCflag == 1)
      {
         yprime[ind1] = yprime[ind1] - vpc.Variable_G_c(a) * y[ind1]; 
         yprime[ind2] = yprime[ind2] - vpc.Variable_G_c(a) * y[ind2];          
         }

      for(int l=0;l<=Variables.lmaxt;l++)
      {
	      yprime[ind1+l]=0.0;
         yprime[ind2+l]=0.0;
         if(Variables.VPCflag == 1)
         {
            yprime[ind1+l] = yprime[ind1+l] - vpc.Variable_G_c(a) * y[ind1+l]; 
            yprime[ind2+l] = yprime[ind2+l] - vpc.Variable_G_c(a) * y[ind2+l];          
            }         
         }
 		}

	/*********************************************************************************************************************************/
	/*  Massless Neutrino perturbations																																*/
	/*********************************************************************************************************************************/
	yprime[ind3]=-Variables.ak*y[ind3+1]-htpr;
   
   if(Variables.VPCflag == 1)
      yprime[ind3] = yprime[ind3] - vpc.Variable_G_c(a) * y[ind3]; 

   for(int l=1;l<=Variables.lmaxnr-1;l++)
   {
   	yprime[ind3+l]=Variables.ak*(l*y[ind3-1+l]-(l+1)*y[ind3+1+l])/(2*l+1);					// Ma, Bertschinger Eq(49)
      if(Variables.VPCflag == 1)
         yprime[ind3+l] = yprime[ind3+l] - vpc.Variable_G_c(a) * y[ind3+l]; 
      }
																											// Truncate moment expansion
	yprime[ind3+Variables.lmaxnr]=Variables.ak*y[ind3-1+Variables.lmaxnr]-(Variables.lmaxnr+1)/tau*y[ind3+Variables.lmaxnr];			// Ma, Bertschinger Eq(49), Eq(51) Combined
   if(Variables.VPCflag == 1)
      yprime[ind3+Variables.lmaxnr] = yprime[ind3+Variables.lmaxnr] - vpc.Variable_G_c(a) * y[ind3+Variables.lmaxnr]; 
   }
   else
   {
      double beta = Variables.ak;
      double betar = beta*Variables.r;

      if (tcp == 1)
      {            
                                                                                                                     // No tight coupling.
           yprime[ind1+2]=beta*(-Curvature.bt1(3,betar)*y[ind1+3]/7.0)-opac*y[ind1+2] + sourcet;                     // Temperature
           yprime[ind2+2]=beta*(-2.0*y[ind3+2]/3.0-Curvature.bt2(3,betar)/7.0*y[ind2+3])-opac*y[ind2+2] + sourcee;   // E polarization
           yprime[ind3+2]=beta*( +2.0*y[ind2+2]/3.0-Curvature.bt2(3,betar)/7.0*y[ind3+3])-opac*y[ind3+2] + sourceb;  // B polarization

          
            if(Variables.VPCflag == 1)
            {
               yprime[ind1+2] = yprime[ind1+2] - vpc.Variable_G_c(a) * y[ind1+2]; 
               yprime[ind2+2] = yprime[ind2+2] - vpc.Variable_G_c(a) * y[ind2+2]; 
               yprime[ind3+2] = yprime[ind3+2] - vpc.Variable_G_c(a) * y[ind3+2];                         
               }   

          for(int l=3;l<=Variables.lmaxt-1;l++)
          {
            yprime[ind1+l]=beta*(Curvature.bt1(l,betar)*y[ind1+l-1]/(2.0*l-1.0) -Curvature.bt1(l+1,betar)*y[ind1+l+1]/(2.0*l+3.0))-opac*y[ind1+l]; // Temperature
            yprime[ind2+l]=beta*(Curvature.bt2(l,betar)/(2.0*l-1.0)*y[ind2+l-1]-4.0/l/(l+1.0)*y[ind3+l]-Curvature.bt2(l+1,betar)/(2.0*l+3.0)*y[ind2+l+1])-opac*y[ind2+l]; // E polarization
            yprime[ind3+l]=beta*(Curvature.bt2(l,betar)/(2.0*l-1.0)*y[ind3+l-1]+4.0/l/(l+1.0)*y[ind2+l]-Curvature.bt2(l+1,betar)/(2.0*l+3.0)*y[ind3+l+1])-opac*y[ind3+l]; // B polarization

            if(Variables.VPCflag == 1)
            {
               yprime[ind1+l] = yprime[ind1+l] - vpc.Variable_G_c(a) * y[ind1+l]; 
               yprime[ind2+l] = yprime[ind2+l] - vpc.Variable_G_c(a) * y[ind2+l]; 
               yprime[ind3+l] = yprime[ind3+l] - vpc.Variable_G_c(a) * y[ind3+l];                         
               }   
          }

         //double xc=tau/r;
         // betar=beta*r;
         // alpha=Curvature.coshK(xc)/Curvature.sinhK(xc);

         // close the hierarchy
         // Just set to zero last multipole. Fancier scheme is not working well.
         // This is faster and still accurate.

         yprime[Variables.lmaxt+ind1]=0.0;
         yprime[ind2+Variables.lmaxt]=0.0;
         yprime[ind3+Variables.lmaxt]=0.0;
         if(Variables.VPCflag == 1)
         {
            yprime[ind1+Variables.lmaxt] = yprime[ind1+Variables.lmaxt] - vpc.Variable_G_c(a) * y[ind1+Variables.lmaxt];
            yprime[ind2+Variables.lmaxt] = yprime[ind2+Variables.lmaxt] - vpc.Variable_G_c(a) * y[ind2+Variables.lmaxt];
            yprime[ind3+Variables.lmaxt] = yprime[ind3+Variables.lmaxt] - vpc.Variable_G_c(a) * y[ind3+Variables.lmaxt];
            }   
         }
      else
      {
         // tight coupling
         shearg=-4.0*htpr/opac/3.0;
         //e2=-sqrt(6.0)/4.0*shearg;
         y[ind1+2]=-4.0*htpr/opac/3.0;
         y[ind2+2]=-sqrt(6.0)/4.0*shearg;
         y[ind3+2]=0.0;
         

         // rest is 0
         yprime[ind1+2]=0.0;
         yprime[ind2+2]=0.0;
         yprime[ind3+2]=0.0;

         if(Variables.VPCflag == 1)
         {
            yprime[ind1+2] = yprime[ind1+2] - vpc.Variable_G_c(a) * y[ind1+2];
            yprime[ind2+2] = yprime[ind2+2] - vpc.Variable_G_c(a) * y[ind2+2];
            yprime[ind3+2] = yprime[ind3+2] - vpc.Variable_G_c(a) * y[ind3+2];
            }   

         for(int l=3;l<=Variables.lmaxt-1;l++)
         {
            yprime[l+ind1]=0.0;
            yprime[ind2+l]=0.0;
            yprime[ind3+l]=0.0;
            if(Variables.VPCflag == 1)
            {
               yprime[ind1+l] = yprime[ind1+l] - vpc.Variable_G_c(a) * y[ind1+l]; 
               yprime[ind2+l] = yprime[ind2+l] - vpc.Variable_G_c(a) * y[ind2+l]; 
               yprime[ind3+l] = yprime[ind3+l] - vpc.Variable_G_c(a) * y[ind3+l];                         
               }   
            }

         }
      }

	}


/************************************************************************************************************************************/
/* This part is based on 																																				*/
/*	1) Signatures of a Graviton Mass in the Cosmic Microwave Background																					*/
/* 	:: Sergei Dubovsky, Raphael Flauger, Alexei Starobinsky, Igor Tkachev																			*/
/*	2)	Microwave Background Constraints on Cosmological Parameters																							*/
/*		:: Matias Zaldarriaga, David N. Spergel, Uro�s Seljak																									*/
/*	3) Signature of Gravity Waves in Polarization of the Microwave Background 																			*/
/*		:: Uro�s Seljak, Matias Zaldarriaga																															*/
/* It will calculate the tensor source term for the CMBR 																									*/
/************************************************************************************************************************************/

void foutputt(int n,double y[],double ypr[],int j,double tau0,double tau,double xsanpass[])
{
	double x,x2,htpr,htdpr,psie,psiedot,psieddot;
	double dt1,dte1,dtb1;
	int ind1,ind2,ind3;

   curvature Curvature(-Variables.omegak);

   if(fabs(Variables.omegak)<0.0001){
      x=Variables.ak*(tau0-tau);
      x2=x*x;
      ind1=4;
      ind2=ind1+Variables.lmaxt+1;
      htpr=y[3];
      htdpr=ypr[3];

                                                                                                // Sergei Eq(8), Zaldarriaga Eq(6)
      psie=y[ind1]/10.0+y[ind1+2]/7.0+3.0*y[ind1+4]/70.0-3.0*y[ind2]/5.0+6.0*y[ind2+2]/7.0-3.0*y[ind2+4]/70.0;
      psiedot=ypr[ind1]/10.0+ypr[ind1+2]/7.0+3.0*ypr[ind1+4]/70.0-3.0*ypr[ind2]/5.0+6.0*ypr[ind2+2]/7.0-3.0*ypr[ind2+4]/70.0;
      psieddot=-0.3*(GlobalArray.opac[j]*psiedot+GlobalArray.dopac[j]*psie)-0.1*htdpr-Variables.ak*(3.0*ypr[ind1+1]/70.0+ypr[ind1+3]/15.0+ypr[ind1+5]/42.0-33.0*ypr[ind2+1]/35.0+8.0*ypr[ind2+3]/15.0-ypr[ind2+5]/42.0);

      if (x>0.0)
      {
         dt1=(-GlobalArray.expmmu[j]*htpr+GlobalArray.vis[j]*psie)/x2;																	// Zaldarriaga Eq(6)
         dte1=GlobalArray.vis[j]*(psie-psieddot/Variables.ak2-6.0*psie/x2-4.0*psiedot/Variables.ak/x)-GlobalArray.dvis[j]*(4.0*psie/x/Variables.ak+2.0*psiedot/Variables.ak2)-GlobalArray.ddvis[j]*psie/Variables.ak2;
         dtb1=2.0*(GlobalArray.vis[j]*(2.0*psie/x+psiedot/Variables.ak)+GlobalArray.dvis[j]*psie/Variables.ak);										// Zaldarriaga Eq(6)
         }
      else
      {
         dt1=0.0;
         dte1=0.0;
         dtb1=0.0;
         }

      dte1=-dte1;
      dtb1=-dtb1;
   }
   else
   {
      double sqrt6 = sqrt(6.0);
      double chi=(tau0-tau)/Variables.r;
      double sinhchi=Curvature.sinhK(chi);
      double rsinh=Variables.r*sinhchi;
      double coshchi=Curvature.coshK(chi);
      double coshchi2=coshchi*coshchi;
      double rsinh2=rsinh*rsinh;
      double beta = Variables.ak;
      // Tensors
      ind1=2;
      ind2=ind1+Variables.lmaxt-1;
      ind3=ind2+Variables.lmaxt-1;
      htpr=y[3];
      htdpr=ypr[3];

      double betar=beta*Variables.r;
      double beta2 = beta*beta;

      double p2=0.1*(y[ind1+2]-sqrt6*y[ind2+2]);
      double p2dot=0.1*(ypr[ind1+2]-sqrt6*ypr[ind2+2]);

      double sourcetdot=GlobalArray.opac[j]*p2dot+GlobalArray.dopac[j]*p2-htdpr;
      double sourceedot=-(sqrt6*GlobalArray.dopac[j]*p2+sqrt6*GlobalArray.opac[j]*p2dot);

      double t2ddot=beta*(-Curvature.bt1(3,betar)*ypr[ind1+3]/7.0)-GlobalArray.opac[j]*ypr[ind1+2]-GlobalArray.dopac[j]*y[ind1+2] + sourcetdot;
      double e2ddot=beta*( -2.0*ypr[ind3+2]/3.0-Curvature.bt2(3,betar)/7.0*ypr[ind2+3])-GlobalArray.opac[j]*ypr[ind2+2]-GlobalArray.dopac[j]*y[ind2+2]+ sourceedot;

      double p2ddot=0.1*(t2ddot-sqrt6*e2ddot);

//beta*( -2.0*ypr[ind3+2]/3.0
      //printf("\nP2 %d %e %e %e",ind2+2,y[ind1+2],sqrt6,y[ind2+2]);
      //printf("\n%e %e %e %e %e",GlobalArray.opac[j],p2dot,GlobalArray.dopac[j],p2,htdpr);
      //printf("\n%d %e %e %e %e",j,beta,betar, Curvature.bt1(3,betar),sourcetdot);
      //printf("\n\n %e %e",beta*( -2.0*ypr[ind3+2]/3.0-Curvature.bt2(3,betar)/7.0*ypr[ind2+3])-GlobalArray.opac[j]*ypr[ind2+2]-GlobalArray.dopac[j]*y[ind2+2], sourceedot);
      //printf("\n\n%e %e %e",t2ddot,e2ddot,p2ddot); exit(1);
      double abeta=Variables.r*Variables.r/sqrt(betar*betar-Variables.kcurv*4.0)/sqrt(betar*betar-Variables.kcurv);

      x=beta*rsinh;
      if (x > 1.0e-8){
         dt1=abeta*(-GlobalArray.expmmu[j]*htpr+GlobalArray.vis[j]*p2)/rsinh2;
         dte1=(GlobalArray.ddvis[j]*p2+2.0*GlobalArray.dvis[j]*p2dot+GlobalArray.vis[j]*p2ddot) +4.0*coshchi*(GlobalArray.dvis[j]*p2+GlobalArray.vis[j]*p2dot)/rsinh + ( -beta2+3.0/Variables.r/Variables.r+6.0*coshchi2/rsinh2)*GlobalArray.vis[j]*p2;
         dte1=dte1*abeta;
         dtb1=(GlobalArray.dvis[j]*p2+GlobalArray.vis[j]*p2dot+2.0*coshchi/rsinh*GlobalArray.vis[j]*p2);
         dtb1=2.0*beta*abeta*dtb1;
         }
      else
      {
         dt1=0.0;
         dte1=0.0;
         dtb1=0.0;
         }  
   }
  // printf("\n%e %e %e",dt1,dte1,dtb1);
   xsanpass[0] = dt1;
   xsanpass[1] = dte1;
   xsanpass[2] = dtb1;
   }

/*************************************************************************************************************/
/* 	This subroutine computes the power spectra for mode ak of the tensor perturbations.							 */
/*************************************************************************************************************/

void powertflat(double ak,int in,double *apower)
{
	double aknlog;
	aknlog=log(ak/0.05);
   *apower=exp((Variables.ant[in]+.5*Variables.alphant[in]*aknlog)*aknlog);
	}

void powertopen(double ak,int in,double *apower)
{
   // This subroutine computes the power spectra
   // for mode ak of the scalar perturbations in the open case. 
   // Now it is set to a power law in k which is slightly more
   // complicated for beta(=ak).

   // Normalize so that tilt leaves the smae power 

   double anorm=0.002;
   double anorm2=anorm*anorm;
   double akn=ak/anorm;
   double pi=3.14159265;
   double aux=tanh(pi*ak*Variables.r/2.0);
   *apower = (ak*ak*Variables.r*Variables.r-4.0*Variables.kcurv)/(-Variables.kcurv+ak*Variables.r*ak*Variables.r)*aux/akn*exp(Variables.ant[in]*0.5*log((-3.0*Variables.kcurv/Variables.r/Variables.r+ak*ak)/anorm2));
   *apower = *apower*2.0e4;
   }


};

void Fderivst(int n,double x,double y[],double yprime[])
{
    tensor Tensor;
    Tensor.fderivst(n,x,y,yprime);
}

#endif