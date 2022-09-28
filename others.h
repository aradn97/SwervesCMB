#ifndef	_OTHERS_H_
#define	_OTHERS_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "variables.h"
#include "initval.h"
#include "darkenergy.h"
#include "neutrino.h"
#include "recombination.h"
#include "doubleinflation.h"

/********************************************************************************************************************/
/*	  This funtion will calculate the inverse of the a_dot. We can intigrate this function with respect to a to      */
/*   get the conformal time at some point																									  */	
/*	  This fuction will use the function nu1(),dynrho() for calculating neutrino density									  */
/********************************************************************************************************************/

double dtauda(double a)
{
	double rx[2];
	double rhonu;
	double dtauda1;
	double grho2; 																								// Its a tempurary veriable for storing the total density
	
	if (a<1.0e-9)
   {
	   rhonu=1.0;
      }
   else
   {
      nu1(a,rx);																								// nu1() will calculate the density of massive neutrinos
      rhonu=rx[0];																							// by interpolating the values already calculated.
      }

																													// 8*pi*G*rho*a**4.
   if (ndyn==0)
	   omegavdyn = omegav;
   else
      omegavdyn = omegav*dynrho(a);


   	grho2=grhom*(omegac+omegab)*a+grhog+grhor*annur+grhonr*annunr*rhonu+grhom*omegavdyn*pow(a,4);//+grhom*omegak*a*a;

	if(dimflag==1)
   	grho2=pow((sqrt(grho2)+sqrt(omegav*grhom*pow(a,4))),2)+grhom*omegak*a*a;

   dtauda1=sqrt(3.0/grho2);
   return(dtauda1);
   }	


/************************************************************************************************************************************/
/*		It will calculate derivative of sound horizon with respect to a.It can be intigated to find the sound horizon.						*/
/************************************************************************************************************************************/

double dsoundda(double a)
{
	double cs,R;
	double dtauda1,dsoundda1;
	
	/*********************************************************************************************************************************/
	/*   Compute expansion rate.																																		*/
	/*********************************************************************************************************************************/

	dtauda1=dtauda(a);
	
	//R = 3*grhob*a / (4*grhog)
   
	R=30000*a*omegab*pow((h0/100.0),2); 
   cs=1.0/sqrt(3.0*(1.0+R));
   dsoundda1=dtauda1*cs;
	return(dsoundda1);
	}	

/************************************************************************************************************************************/
/*		finital will calculate the initial values for the source terms depending on the choosen initial condition							*/
/*		It will use the functions : dynrho(),																														*/
/************************************************************************************************************************************/

void finitial(double y[],double tau)
{
	double rhosan[2];
	double rhonu,pnu;
	double a,a2;
	
	double grho,fracnu,s,gpres;
	double psi,C,akt2,h,deltac,deltag,deltab,deltar;
	double thetac,thetag,thetab,thetar,shearr,ahdot;
	double dq,q,aq,v,akv,f1,yrad,weos;
	double rc_phi,rc_psi,expq,eta;
	double deltan,thetan,delta0;
	double dlfdlq1;
	
	int ind;
	
	a=tau*adotrad;																// a at the initial point	
   a2=a*a;	
							
   nu1(a,rhosan);																// pressure and density at initial point   
   rhonu=rhosan[0];
   pnu=rhosan[1];
   
   omegavdyn = omegav*dynrho(a);
   weos = wdyn_func(a);
																					
	if(dimflag==0)																// 8*pi*G*rho*a**2. energy density for FRW and 5 dimensions
   	grho=grhom*(omegac+omegab)/a+(grhog+grhor*annur+grhonr*annunr*rhonu)/a2+grhom*omegavdyn*a2;
	else
	{
		grho=sqrt(grho)+sqrt(omegav*grhom*a2);
		grho=grho*grho;
		}
																					// pressure 8*pi*G*P*a**2.
	gpres=((grhog+grhor*annur)/3.0+grhonr*annunr*pnu)/a2+weos*grhom*omegavdyn*a2;
																					
   s=grho+gpres;
   fracnu=(4.0/3.0)*(grhor*annur+grhonr*annunr)/(a2*s);			// Nutrino Fraction

																					// Use yrad=rho_matter/rho_rad to correct initial conditions for matter+radiation.  
   yrad=grhom*(omegac+omegab)*a/(grhog+grhor*annur+grhonr*annunr*rhonu);
   
   
	/**********************************************************************************************************/   
	/* Choose one of the following four cases for initial conditions, or add your own.                        */  
	/**********************************************************************************************************/   
   
	/***********************************************************************************************************/   
	/*  First case.																														  */
	/*  Isentropic ("adiabatic") initial conditions. normalize to zeta=1													  */
	/*  Check the equation set (96),(97),(98) Ma,Bertschinger                           							  */
	/***********************************************************************************************************/	
   
   if (initfl==1)
   {
	   psi=-1.0/(1.5+2.0*fracnu/5.0);													 //fracnu = rhonu/(rhogama+rhonu)	Eq(98f)				
		C=(15.0+4.0*fracnu)/20.0*psi;                                           
      akt2=(ak*tau)*(ak*tau);																  	
      h=C*akt2*(1.0-0.2*yrad);															 //Eq(96g)	
      eta=2.0*C-(5.0+4.0*fracnu)*C*akt2*(1.0-yrad/3.0)/(6.0*(15.0+4.0*fracnu));   //Eq(96h)   
      f1=(23.0+4.0*fracnu)/(15.0+4.0*fracnu);								       //Eq(96e f1=thetanu/thetagamma)
      deltac=-0.5*h;																			 //Eq(96a, 96b, 96g)
      deltag=-2.0*h*(1.0-akt2/36.0)/3.0;  											 //Eq(96a)
      deltab=0.75*deltag;																	 //Eq(96b)
      deltar=-2.0*h*(1.0-akt2*f1/36.0)/3.0;   										 //Eq(96a, 96b)
      thetac=0.0;																				 //Eq(96c)											
      thetag=-C*akt2*akt2/(tau*18.0);													 //Eq(96d)	
      thetab=thetag;																			 //Eq(96d)			
      thetar=f1*thetag;																		 //Eq(96e)	
      shearr=(4.0/15.0)*(ak2/s)*psi*(1.0+7.0/36.0*yrad); 						 //Eq(98d)
      ahdot=2.0*C*ak2*tau*a*(1.0-0.3*yrad);											 //Eq(94a) a*hdot
		}
		
	/***********************************************************************************************************/
	/*  Second case.																														  */
	/*  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.												  */
	/*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois									  */
	/*  Check Eq (7),(8),(9)     Sb=0    Snu=0																					  */
	/***********************************************************************************************************/
   else if(initfl==2)																					
   {																									
	   delta0=1.0;																									
	   h=delta0*yrad*(1.0/(1.0+omegab/omegac)-0.5*yrad);															
      deltac=delta0-0.5*h;																	// Perturb CDM	
      deltag=-2.0/3.0*h;																	//
      deltab=0.75*deltag;																	// Sb = 0
      deltar=deltag;																			// Snu = 0	
      thetac=0.0;
      thetag=-h*ak2*tau/12.0;
      thetab=thetag;																			//Eq(11, 7)
      thetar=thetag;																			//Eq(11)
      shearr=0.0;
      ahdot=adotrad*h*(1.0-0.5*yrad);
      eta=-h/6.0;
      }

	/***********************************************************************************************************/
	/*  Third case.																														  */
	/*  Isocurvature CDM initial conditions: perturb only baryons as a --> 0.											 */
	/*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois									  */
	/*  Check Eq (7),(8),(9)     Sr=0    Snu=0					   															  */
	/***********************************************************************************************************/      
   else if(initfl==3)
   {
	   delta0=1.0;
      h=delta0*yrad*(1.0/(1.0+omegac/omegab)-0.5*yrad);
      deltab=delta0-0.5*h;
      deltac=-0.5*h;
      deltag=-2.0/3.0*h;
      deltar=deltag;
      thetac=0.0;
      thetag=-h/12.0*ak2*tau;
      thetab=thetag;
      thetar=thetag;
      shearr=0.0;
      ahdot=adotrad*h*(1.0-0.5*yrad);
      eta=-h/6.0;
      }
	/***********************************************************************************************************/
	/*  Third case.																														  */
	/*  Isocurvature CDM initial conditions: unperturbed everything 													     */
	/*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois									  */
	/*  Check Eq (7),(8),(9)     Sb=0    Snu=0	Sr =9									    									  */
	/***********************************************************************************************************/   
	else if (initfl==4)
	{
		delta0=1.0;
      h=delta0*yrad*(1.0/(1.0+omegac/omegab)-0.5*yrad);
      deltab=-0.5*h;
      deltac=-0.5*h;
      deltag=-2.0/3.0*h;
      deltar=deltag;
      thetac=0.0;
      thetag=-h/12.0*ak2*tau;
      thetab=thetag;
      thetar=thetag;
      shearr=0.0;
      ahdot=adotrad*h*(1.0-0.5*yrad);
      eta=-h/6.0;
      }

	/***********************************************************************************************************/
	/*	 If some other value then neglect and exit form the program															  */
	/***********************************************************************************************************/
   else
   {
	   printf("Its not possible. Initial condition is outside the range of values");
      exit(1);
      }
	
   deltan=deltar;
   thetan=thetar;

   y[1]=a;
   y[2]=ahdot;
   y[3]=eta;
														  									//CDM.
   y[4]=deltac;
   y[5]=thetac;
																							//Baryons.
   y[6]=deltab;
   y[7]=thetab;
																							//Photons (total intensity and polarization).
   y[8]=deltag;
   y[9]=thetag;
 
  for(int lsa=1;lsa<=lmaxg;lsa++)
   {
	   y[9+lsa]=0.0;																		// Ma Bertschingger Eq(64) F_gamma
      y[9+lmaxg+lsa]=0.0;															   // Ma Bertschingger Eq(64) G_gamma
      }
																							// Massless neutrinos.
   y[10+2*lmaxg]=deltar;																
   y[11+2*lmaxg]=thetar;
   y[12+2*lmaxg]=shearr*2.0;														// Ma Bertschingger Eq(49) F_nu_2
   
   for(int l=3;l<=lmaxnr;l++)
   	y[10+2*lmaxg+l]=0.0;															// Ma Bertschingger Eq(49) F_nu_l		
      
     																						//  Massive neutrinos.
	dq=1.0;
   for(int i=1;i<=nqmax;i++)
   {
		q=i*dq-0.5;
      aq=a*amnu/q;
      v=1.0/sqrt(1.0+aq*aq);
      akv=ak*v;
      expq=exp(-q);
      dlfdlq1=-q/(1.0+expq);
      y[iq0+i-1]=-0.25*dlfdlq1*deltan;											// Ma Bertschingger Eq(97) Psi_nu_0 for different velovities of neutrino (massive)
      y[iq1+i-1]=-dlfdlq1*thetan/(3.0*akv);									// Ma Bertschingger Eq(97) Psi_nu_1 for different velovities of neutrino (massive)
      y[iq2+i-1]=-0.5*dlfdlq1*shearr;											// Ma Bertschingger Eq(97) Psi_nu_3 for different velovities of neutrino (massive)
      for(int l=3;l<=lmaxnu;l++)
      {
	      ind=10+2*lmaxg+lmaxnr+i+l*nqmax;
         y[ind]=0.0;																	// Ma Bertschingger Eq(97) (56) Psi_nu_l for diff velovities of neutrino (massive)
         }
      }
																							//  Check energy constraint equation.
													
   if(ndyn==1 || ndyn==2)														// Perturbations in dark energy
   {
	   rc_phi = 0.0;
      rc_psi = 0.0;
      y[nvar0-1] = rc_phi;
      y[nvar0]   = rc_psi;
		}	
	}


/************************************************************************************************************************************/  
/* 	This will calculate the derivatives of the scaler perturbation terms. This part is totally based on Ma & Bertschunger			*/
/************************************************************************************************************************************/  		
		
void fderivs(int n,double x,double y[],double yprime[])
{																											// n store total number of values in y
	double ep0 = .01;
	double xsanpass[3],ep,tau,epstc,dpnu,grho,gpres;
	double deltac,thetac,deltab,thetab,deltag,thetag,shearg;
	double polter,deltar,thetar,shearr,a2,cs2,dlxedla;										//,opec
	double a,eta,rc_phi,rc_psi,photbar,pb43,ysanpass[2];
	double tcp,tcp1,tcp2,tcp3,drag,deltabdot;
	double rgrho,slip,pnu,fnu,shearnu;
	double dgrho,dgpres,dgtheta,etadot;
	double deltacdot,thetacdot,deltagdot,thetabdot;
	double rhonudot,shearnudot,adotdota,weos;       
   double thetagdot,deltardot,thetardot,Fnu2dot;
   double dgprs_phi,dahdotdtau,sgrhooa2;
   double opac1,dgrho_phi,dgtheta_phi,dq,q,v,aq;
   double akv[nqmax0+1];  //,rc_dphi,rc_dpsi;
   double zsanpass[4];

   int ind,tcpa,tcpb;	
		
	if (ak > epsw) 	 																				// ep is used to stop the tight coupling approximation.
   	ep=ep0;
   else
   	ep=0.5*ep0;																							

	tau=x;																								// time of the derivative
   
   a=y[1];
   eta=y[3];
        
   deltac=y[4];																						//  CDM.
   thetac=y[5];
        
   deltab=y[6];																						// Baryons.
   thetab=y[7];
      
   deltag=y[8];																						//  Photons.
   thetag=y[9];
   shearg=y[10]/2.0;
	
	
	polter=y[10]+y[9+lmaxg]+y[11+lmaxg];															//  Polarization term. Ma Eq(60). F2+G0+G2

   deltar=y[10+2*lmaxg];																				//  Massless neutrinos.
   thetar=y[11+2*lmaxg];
   shearr=y[12+2*lmaxg]/2.0;

   a2=a*a;
   thermo(tau,xsanpass);																			// Get the Value of dlxedla
	
	cs2=xsanpass[0];
	opac1=xsanpass[1];
	dlxedla=xsanpass[2];
																								
   photbar=grhog/(grhom*omegab*a);																// Photon mass density over baryon mass density.
   pb43=(4.0/3.0)*photbar;

	tcp=0.0;																								// Tight Coupling parameters
   tcp1=ak/opac1;													
   tcp2=1.0/(opac1*tau);													
   tcp3=ak*tau;

   tcpa=0;
   tcpb=0;
   
   epstc=max(tcp1,tcp2);
	
   if ((epstc > 5.0e-3) && (a > 1.35e-5))
   	tcpa=1;
   
   epstc=epstc/(1.0+pb43);
   
   if ((epstc > 5.0e-3) && (a > 1.35e-4)) 
   	tcpb=1;

	/*********************************************************************************************************************************/  
	/*  Compute expansion rate (a dot = sqrt(8 pi G rho / 3)).																								*/
	/*********************************************************************************************************************************/ 

   if (amnu==0.0)																						// No massive Neutrino, so no perturbation  			
   {
	   rhonu=1.0;
      pnu=1.0/3.0;
      drhonu=0.0;
      fnu=0.0;
      dpnu=0.0;
      shearnu=0.0;
      }			
   else																									// If massive Neutrino present then perturbe		
   {
	   nu1(a,ysanpass);
	   rhonu=ysanpass[0];
	   pnu=ysanpass[1];
      nu2(a,zsanpass,y+iq0,y+iq1,y+iq2);
      drhonu = zsanpass[0];
      fnu = zsanpass[1];
      dpnu = zsanpass[2];
      shearnu = zsanpass[3];
      }	
	
   omegavdyn = omegav * dynrho(a);
   weos = wdyn_func(a);
	
																								//  8*pi*G*rho*a^2 and 8*pi*G*P*a^2.
   
   grho=grhom*(omegac+omegab)/a+(grhog+grhor*annur+annunr*rhonu*grhonr)/a2+grhom*omegavdyn*a2;
		  																											
	if(dimflag==1)	
   	grho=pow((sqrt(grho)+sqrt(omegav*grhom*a2)),2); 									// energy density for  5 dimensions
																											
	adotoa=sqrt(grho/3.0);																			
   yprime[1]=adotoa*a;																						
																												
	gpres=((grhog+grhor*annur)/3.0+grhonr*annunr*pnu)/a2;									// 8*pi*G*rho*a^2 and 8*pi*G*P*a^2.
   gpres=gpres + weos*grhom*omegavdyn*a2;
   
	if (ndyn==1 || ndyn==2)																			// dark energy perturbations
	{
		rc_phi = y[n-1];
      rc_psi = y[n]; 

      dyn_nrg(a,grho,gpres,rc_phi,rc_psi,xsanpass);
      dgrho_phi = xsanpass[0];
      dgprs_phi = xsanpass[1];
      dgtheta_phi=xsanpass[2];
      }
																											// Evaluate metric and massive neutrino perturbations.
																												//  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
   dgrho=grhom*(omegac*deltac+omegab*deltab)/a +(grhog*deltag+grhor*annur*deltar+grhonr*annunr*drhonu)/a2;
   
   if(ndyn==1 || ndyn==2)
   	dgrho=dgrho+dgrho_phi;
   	
	dgpres=(grhog*deltag+grhor*annur*deltar)/(a2*3.0) +grhonr*annunr*dpnu/a2;
	
	if (ndyn==1 || ndyn==2)
		dgpres=dgpres+dgprs_phi;

   if (initfl==4)																						// Add a seed if desired.
   	dgrho=dgrho+grhom/a;
   
   dahdotdtau=-(dgrho+3.0*dgpres)*a;
   yprime[2]=dahdotdtau;
																											// Force energy conservation.
   hdot=(2.0*ak2*eta+dgrho)/adotoa;																// Eq(21a)
																											// 8*pi*G*(rho+P)*theta*a^2.
   dgtheta=grhom*(omegac*thetac+omegab*thetab)/a+(4.0/3.0)*(grhog*thetag+annur*grhor*thetar)/a2+annunr*grhonr*ak*fnu/a2;
   
   
   if(ndyn==1 || ndyn==2) 
   	dgtheta=dgtheta+dgtheta_phi;

   etadot=0.5*dgtheta/ak2;																			//Eq(21b)
   yprime[3]=etadot;


	if (tcpa!=1) 																						// First order approximation for photon shear
    	shearg=(1.0/(opac1*9.0))*((8.0/3.0)*thetag+(4.0/3.0)*hdot+8.0*etadot);
      

   dgshear=(4.0/3.0)*(grhog*shearg+annur*grhor*shearr)/a2+annunr*grhonr*shearnu/a2;// 8*pi*G*(rho+P)*sigma*a**2.
	

	/*********************************************************************************************************************************/     
	/*	  CDM equations of motion.																																		*/			
	/*********************************************************************************************************************************/
   
	deltacdot=-thetac-0.5*hdot;																	// Ma & Bertschunger Eq. (29), (42) As not ineracting w=0 and delP=0
   yprime[4]=deltacdot;
   thetacdot=-adotoa*thetac;																		// Ma & Bertschunger Eq(29) As not ineracting w=0 and delP=0 K=0 as flat
   yprime[5]=thetacdot;
       
	/*********************************************************************************************************************************/     
	/*  Photon equations of motion (total intensity).																											*/			
	/*********************************************************************************************************************************/
   deltagdot=-(4.0/3.0)*(thetag+0.5*hdot);													// Ma & Bertschunger Eq(63)
   yprime[8]=deltagdot;																				// As deltag in involved in baryon calculation so it is calculated here

	/*********************************************************************************************************************************/     
	/*	  Baryon equations of motion.																																	*/			
	/*********************************************************************************************************************************/     
   deltabdot=-thetab-0.5*hdot;																	// Ma & Bertschunger Eq(66)
   yprime[6]=deltabdot;
     
   drag=opac1*(thetag-thetab);
     
   if(tcpb==1)																							// If baryons and photons are uncoupled.
		thetabdot=-adotoa*thetab+ak2*cs2*deltab+pb43*drag;									// Ma & Bertschunger Eq(66) 
		  	
   else																									//If baryons and photons are tightly coupled. Zeroth-order approximation to baryon velocity.
	{
		thetabdot=(-adotoa*thetab+ak2*cs2*deltab+ak2*pb43*(0.25*deltag-shearg))/(1.0+pb43); //

		if(dimflag==1)
		{
			if (amnu==0.0)
			{
				rhonudot=0.0;
           	shearnudot=0.0;
				}        		
         else
         {
	         nuder(a,rhonu,ysanpass,adotoa,y+iq2-1,yprime+iq2-1);
            rhonudot = ysanpass[0];
            shearnudot = ysanpass[1];
            }
      	
      	sgrhooa2=sqrt(grho/a2);
      	
      	rgrho=sqrt(grhom*omegav*grhom*omegav+grhom*(omegab+omegac)/(a2*a)+(grhog+grhor*(annur+annunr*rhonu))/(a2*a2));
      	adotdota=sgrhooa2/rgrho*(grhom*(omegab+omegac)*(-3/(2*a))+(grhog+grhor*(annur+annunr*rhonu))*(-2/a2)+grhor*annunr*rhonudot/(2*a2*adotoa))+2*grho;
      	
      	adotdota=adotdota/3.0;
      	}
		
		else
		{
			adotdota=0.5*(adotoa*adotoa-gpres);
			}
		
																											//  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
      slip=(2.0*pb43/(1.0+pb43)+dlxedla)*adotoa*(thetab-thetag)+(1.0/opac1)*(-adotdota*thetab-adotoa*ak2*0.5*deltag+ak2*(cs2*deltabdot-0.25*deltagdot))/(1.0+pb43);
																											//  First-order approximation to baryon velocity.
      thetabdot=thetabdot+pb43/(1.0+pb43)*slip;
      }
		
   yprime[7]=thetabdot;

	/*********************************************************************************************************************************/     
	/*	  Photon equations of motion (polarization).  																											*/			
	/*********************************************************************************************************************************/
	thetagdot=(-thetabdot-adotoa*thetab+ak2*cs2*deltab)/pb43 + ak2*(0.25*deltag-shearg);
   yprime[9]=thetagdot;                                        						// Ma & Bertschunger Eq(75)                      
 
   
	if (tcpa==1)																						// Treat baryons and photons as uncoupled.
   {
	   yprime[10]=(8.0/15.0)*thetag-0.6*ak*y[11]-opac1*y[10]+(4.0/15.0)*hdot+(8.0/5.0)*etadot+0.1*opac1*polter;
																											// Polarization equations for l = 0, 1, 2.
    	
      yprime[9+lmaxg]=-ak*y[10+lmaxg]-opac1*y[9+lmaxg]+0.5*opac1*polter;
      yprime[10+lmaxg]=(ak/3.0)*(y[9+lmaxg]-2.0*y[11+lmaxg])-opac1*y[10+lmaxg];
      yprime[11+lmaxg]=ak*(0.4*y[10+lmaxg]-0.6*y[12+lmaxg])-opac1*y[11+lmaxg]+0.1*opac1*polter;
      for(int l=3;l<=lmaxg-1;l++)
      {
	      yprime[8+l]=ak*(1.0/(2*l+1.0))*(l*y[7+l]-(l+1)*y[9+l])-opac1*y[8+l];
         yprime[9+lmaxg+l]=ak*(1.0/(2*l+1.0))*(l*y[8+lmaxg+l]-(l+1)*y[10+lmaxg+l])-opac1*y[9+lmaxg+l];
			}
		}	
	else																									// Treat baryons and photons as tightly coupled (with no polarization).
   {
	   yprime[10]=0.0;
      yprime[9+lmaxg]=0.0;
      yprime[10+lmaxg]=0.0;
      yprime[11+lmaxg]=0.0;
      for(int l=3;l<=lmaxg-1;l++)
      {
      	yprime[8+l]=0.0;
         yprime[9+lmaxg+l]=0.0;
	      }
      }																									// Truncate moment expansion.

   yprime[8+lmaxg]=ak*y[7+lmaxg]-(lmaxg+1)/tau*y[8+lmaxg]-opac1*y[8+lmaxg];
   yprime[9+2*lmaxg]=ak*y[8+2*lmaxg]-(lmaxg+1)/tau*y[9+2*lmaxg]-opac1*y[9+2*lmaxg];
        
	/*********************************************************************************************************************************/     
	/* 	 Massless Neutrino equations of motion.																												*/			
	/*********************************************************************************************************************************/
   deltardot=-(4.0/3.0)*(thetar+0.5*hdot);													// Ma & Bertschunger Eq(49a)  	
   yprime[10+2*lmaxg]=deltardot;
        	
   thetardot=ak2*(0.25*deltar-shearr);															// Ma & Bertschunger Eq(49b) 		
   yprime[11+2*lmaxg]=thetardot;
        
   Fnu2dot=8.0/15.0*thetar-0.6*ak*y[13+2*lmaxg]+4.0/15.0*hdot+8.0/5.0*etadot;
   yprime[12+2*lmaxg]=Fnu2dot;																		// Ma & Bertschunger Eq(49c) 
        
   for(int l=3;l<=lmaxnr-1;l++)
   	yprime[10+2*lmaxg+l]=ak*(l*y[9+2*lmaxg+l]-(l+1)*y[11+2*lmaxg+l])/(2*l+1);   	// Ma & Bertschunger Eq(49d) 
	
																											// Truncate moment expansion. Eq(51)
   yprime[10+2*lmaxg+lmaxnr]=ak*y[9+2*lmaxg+lmaxnr]-(lmaxnr+1)/tau*y[10+2*lmaxg+lmaxnr];

	/*********************************************************************************************************************************/     
	/*	  Equations for dark energy perturbations.																												*/			
	/*********************************************************************************************************************************/
   if((ndyn==1)||(ndyn==2))
   {
	   dyn_phi(a,hdot,grho,gpres,rc_phi,rc_psi,ysanpass);									// Dark Energy perturbation         
      yprime[n-1] = ysanpass[0];																	//rc_dphi;
      yprime[n]   = ysanpass[1];																	//rc_dpsi;
		}


   if(nqmax==0) 
   	return;																				//Return if no massive Neutrinos are present

	/*********************************************************************************************************************************/     
	/*  Massive Neutrino equations of motion.																														*/			
	/*********************************************************************************************************************************/

	dq=1.0;
	
	for(int i=1;i<=nqmax;i++)
   {
	   q=i*dq-0.5;
      aq=a*amnu/q;
      v=1.0/sqrt(1.0+aq*aq);
      akv[i]=ak*v;
      }	
																											
   for(int i=1;i<=nqmax;i++)
   {
	   ind=iq0+i-1;
      yprime[ind]=-akv[i]*y[ind+nqmax]+hdot*dlfdlq[i]/6.0;								// Ma & Bertschunger Eq(56a) 						
      ind=iq1+i-1;
      yprime[ind]=akv[i]*(y[ind-nqmax]-2*y[ind+nqmax])/3;								// Ma & Bertschunger Eq(56b)
      ind=iq2+i-1;																					// Ma & Bertschunger Eq(56c)			
      yprime[ind]=akv[i]*(2*y[ind-nqmax]-3*y[ind+nqmax])/5 -(hdot/15.0+2.0/5.0*etadot)*dlfdlq[i];
      ind=10+2*lmaxg+lmaxnr+i+lmaxnu*nqmax;													
      																									// Truncate moment expansion.
      yprime[ind]=akv[i]*y[ind-nqmax]-(lmaxnu+1)/tau*y[ind];							// Ma & Bertschunger Eq(57c) Eq(58) Combine
		}	

   for(int l=3;l<=lmaxnu-1;l++)
   {
   	for(int i=1;i<=nqmax;i++)
      {
	      ind=10+2*lmaxg+lmaxnr+i+l*nqmax;
         yprime[ind]=akv[i]*(1.0/(2*l+1.0))*(l*y[ind-nqmax]-(l+1)*y[ind+nqmax]);	// Ma & Bertschunger Eq(57c)
			}
		}
	}	

/*************************************************************************************************************/
/*		This subroutine computes the power spectra for mode ak of the scalar perturbations.							 */		
/*************************************************************************************************************/

void powersflat(double ak,int in,double *apower)
{
	double win,aknlog;
	if (kcutflag==0)
		win=1.0;
   else
   {
	   win=2.0*pow((ak/aksplit),4);
      win=2.0*exp(-win)/(1.0+exp(-win));
      if (kcutflag==-1) win=1.0-win;
      }
	
// Normalize so tilt does not change power at pivot point k=0.05/Mpc
	aknlog=log(ak/0.05);
   *apower=exp((an[in]-1.0+.5*alphans[in]*aknlog + dalphansdlnk[in]*aknlog*aknlog/6.)*aknlog);
//	printf("in \t%e\t%e\t%e\t%e\t%e   out",aknlog,an[in],alphans[in],dalphansdlnk[in],*apower);
// K split
   *apower=(*apower)*win;
//	if(twofield == 1)	
//		*apower = genphiadiabatic(ak);

	}


#endif