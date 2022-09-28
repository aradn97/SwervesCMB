/********************************************************************************************************************/
/*									Cosmic Microwave Background Anisotropy - Numerical Simulation	(CMBANS) 					  */
/********************************************************************************************************************/
/*                                                                                                                  */
/*                                                                                                                  */
/*      The program can calculate the CMB anisotropy and polarization based on the STANDARD cosmological model.     */
/*      The program is based on the CMBFAST code. The code will integrate the linearised perturbation equations     */
/*      numericaly and calculate the theoratical power spectrum of CMB.                                             */
/*                                                                                                                  */
/*      The code is written by Santanu Das during his graduate school project in IUCAA.                             */
/*                                                                                                                  */
/*                                                                                                                  */
/********************************************************************************************************************/


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
#include "recfast.h"
#include "lensing.h"

#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0
#define MAX 50											// For dverk

#define lmaxt0 10										//Tensor equations will be discritize in these many steps
#define lmx0 30										//		
#define nk0 300										//Sets the dimension of the k arrays where the sources will be calculated 

#define nvar0t 2*(lmaxt0+1)+2+1+lmaxnr0+1+1
#define nq0 2000										//Maximum number of element in the (a,w) table for dark energy
#define nkmax 2*l0max
#define ketamax0 3*l0max+126		
#define lm1 8
#define lm2 7
#define ntau 201					 // Sets the sampling time for lensing
/**********************************************************************************************************************************/
/*		All the external stryctures are defined here.																										 */
/**********************************************************************************************************************************/


int l[lmax+1],l0;								// All the l values											/lvalues1/

int lmoin;										// maximum l values to be calculated		

double ztf[ntfmax+1];						//																	/transfer/
int ict,nlnkt;									//	ict = flag CMB/TF : ntf = number of TF				/transfer/

int nflag_rho;									//																	/qstore/

double ant[nnmax+1],rat[nnmax+1];				//																	/tensor/
double alphant[nnmax+1];							//	Flag scaler tensor										/tensor/
int lmaxt,itflag;
int rcrat;												//																	/qtens/

double denl[lmx0+1];			//	denl=1/(2*j+1)												/store/															   
															   
int nvar;												// dark matter possition in the array					/dem/

int nvart;												// Tensor discritisation array							/tensor/
double aa,amnu1;										//																	/nu1d/
double phik[nstep0+1][ntau+1];				// Tau0 will store the present conformal time		/par/
double tol=1.0e-8;									// The tolerance value for intigrations.
int riflag;												// reionization flag (y/n) rc flag (peebles/recfast)
double d[nk0+1][nstep0+1];							//	Source terms 												/memory1/	
double dpr[nk0+1][nstep0+1];						//	Scaler														/memory1/
double dp[nk0+1][nstep0+1];						//																	/memory1/
double dppr[nk0+1][nstep0+1];						//																	/memory1/
double dkk[nk0+1][nstep0+1];						//																	/memory1/	
double dkpr[nk0+1][nstep0+1]; 					//																	/memory1/

double dkpradb[nk0+1][nstep0+1]; 				//																	/memory4/
double dkpriso[nk0+1][nstep0+1]; 				//																	/memory4/

double dadb[nk0+1][nstep0+1];						//	Adiabatic perturbation									/memory3/
double dpadb[nk0+1][nstep0+1];					//																	/memory3/
double dkkadb[nk0+1][nstep0+1];					//																	/memory3/
double diso[nk0+1][nstep0+1];						//	Isocurvature perturbation								/memory3/
double dpiso[nk0+1][nstep0+1];					//																	/memory3/
double dkkiso[nk0+1][nstep0+1];					//																	/memory3/

double dpradb[nk0+1][nstep0+1];					//	Adiabatic perturbation									/memory3/
double dppradb[nk0+1][nstep0+1];					//																	/memory3/
double dkkpradb[nk0+1][nstep0+1];				//																	/memory3/
double dpriso[nk0+1][nstep0+1];					//	Isocurvature perturbation								/memory3/
double dppriso[nk0+1][nstep0+1];					//																	/memory3/
double dkkpriso[nk0+1][nstep0+1];				//																	/memory3/

double dt[nk0+1][nstep0+1];						//	Source terms												/memory2/
double dtpr[nk0+1][nstep0+1];						//	Tensor														/memory2/
double dte[nk0+1][nstep0+1];						//																	/memory2/	
double dtepr[nk0+1][nstep0+1];					//																	/memory2/
double dtb[nk0+1][nstep0+1];						//																	/memory2/
double dtbpr[nk0+1][nstep0+1];					//																	/memory2/
																																
double aw[nq0+1],wd[nq0+1],dwda[nq0+1];		//																	
double wdpr[nq0+1],logrho[nq0+1]; 				//																	/wvar/
   
double power[nk0*6+1][ntfmax+1];					//	Store the transfer function							/transfer_function/ 
int istore[ntfmax+1];								//																	/transfer_function/

double ctl[lmax+1][nnmax+1];						//	
double ctlpr[lmax+1][nnmax+1];					//	Tensor Cl's
double ctel[lmax+1][nnmax+1];						//	
double ctelpr[lmax+1][nnmax+1];					//
double ctbl[lmax+1][nnmax+1];						//
double ctblpr[lmax+1][nnmax+1];					//
double ctcl[lmax+1][nnmax+1];						//
double ctclpr[lmax+1][nnmax+1];					//

int kmaxjl;
double ajl[ketamax0+1][lmax+1];					//
double ajlpr[ketamax0+1][lmax+1],xx[ketamax0+1];	//	BESSEL'S Function
double dxx[ketamax0+1];										//    																	/jlgen/								   

double sand[nk0+1],sandp[nk0+1];					 // tempurary variables for passing	
double sandpr[nk0+1],sandppr[nk0+1];			 // in the splder function	
double sandkpr[nk0+1],sandkk[nk0+1];			 //
        		
double sandt[nk0+1],sandte[nk0+1];				 // tempurary variables for passing
double sandtb[nk0+1],sandtpr[nk0+1];			 // in the splder function	
double sandtepr[nk0+1],sandtbpr[nk0+1];	    //

double bpcl[lmax+1],bpcvl[lmax+1];				 //	
double bpccl[lmax+1],bpckkl[lmax+1];			 // tempurary variables for passing the array
double bpctkl[lmax+1],bpclpr[lmax+1];			 //
double bpcvlpr[lmax+1],bpcclpr[lmax+1];		 //
double bpckklpr[lmax+1],bpctklpr[lmax+1];     //
double bpcplpr[lmax+1],bpcpl[lmax+1];			 //
			
double bpctl[lmax+1],bpctel[lmax+1];			 //	
double bpctbl[lmax+1],bpctcl[lmax+1];			 // tempurary variables for passing the array
double bpctlpr[lmax+1],bpctelpr[lmax+1];		 //	
double bpctblpr[lmax+1],bpctclpr[lmax+1]; 	 //				


/************************************************************************************************************************************/
/*		This function will define the l values for which the Cl will be calculated. For the Other values of l Cl will be				 	*/
/*		Calculated by fitting an spline 																																*/
/************************************************************************************************************************************/

void initlval(int lmoin)
{
	lmo=lmoin+300;
	int lind=1;
	int lvar;
   for(lvar=2;lvar<=10;lvar++,lind++)
      l[lind]=lvar;
      
   l[lind]=12;  lind++;																							
   l[lind]=15;  lind++;
   l[lind]=20;  lind++;
   l[lind]=30;  lind++;
   l[lind]=40;  lind++;
   l[lind]=50;  lind++;
   l[lind]=60;  lind++;
   l[lind]=70;  lind++;
   l[lind]=80;  lind++;
   l[lind]=90;  lind++;
   l[lind]=110; lind++;
   l[lind]=130; lind++;

   for(lvar=150;lvar<=lmo;lvar+=50,lind++)
     l[lind]=lvar;
 
	l0=--lind;
	}
	

/**************************************************************************************************************/
/*		Calculate the massive neutrino density and pressure at different time.   										  */
/*		Called from : initnul()																								           */
/*		It will use : splint()																									    	  */
/**************************************************************************************************************/

// (Fn2) : Massive neutrino P and Rho

void ninul(double a,double rhonu[2])
{
	double qmax=30.0;
	int nq=1000;
	int nq1=1001;
	double dum1[nq1+1],dum2[nq1+1];																// Dummey variable for calculation of pnu and rhonu
	double v,q;
	double qdn1;
	
	if(amnu==0.0)
	{
		rhonu[0]=1.0;
		rhonu[1]=1.0/3;
		}
	
	double dq=qmax/nq,aq;
	dum1[1]=0.0;
	dum2[1]=0.0;
	
	for(int i=1;i<=nq;i++)
	{
		q=i*dq;
		aq=a*amnu/q;
		v=1/sqrt(1.0+aq*aq);
		qdn1=dq*q*q*q/(exp(q)+1.0);
      dum1[i+1]=qdn1/v;
      dum2[i+1]=qdn1*v;
		}
	rhonu[0]=splint(dum1,nq1);
	rhonu[1]=splint(dum2,nq1);
	
	rhonu[0]=(rhonu[0]+dum1[nq1]/dq)/constan;
   rhonu[1]=(rhonu[1]+dum2[nq1]/dq)/(constan*3.0);
	
	}
		
		
/********************************************************************************************************************/
/*	  Initiallizing the massive neutrino pressure and density.                                                       */
/*   The function will then calculate : 1) p 2) rho 3) dp 4) drho 5) ddp 6) ddrho 7) qdn                         	  */
/*   This will be called from : 1) camflat() - P3                                                                   */
/*	  This function will use  : 1) ninul()	 2) splini()  3) splder()																  */
/********************************************************************************************************************/

// (Fn 1) : initiallize pressure and density of massive neutrinos	
	
void initnul()
{  
	double rhonu[2],a,amnu1=1.0e9;
	int yy;
	double q,dq;
	int i;
	FILE *ff;
	
	amin=1.0e-9;
	
	printf("Is there any file containing Nutrino density and pressure (1/0): ");
	scanf("%d",&yy);
	//yy=getch();
	printf("%c",yy);
   
	if(yy==1)
	{
		ff=fopen("neutrino.dat","r");
		if(ff==NULL)
			printf("Sorry the file can not be opened");
				
		fscanf(ff,"%lf",&amnu1);																			// This previous simulation may be for some other neutrino mass
																													// So we need to check that before using it. :)
		if(amnu1/amnu-1.0 < 1e-8)
		{																															
			for(i=1;i<=nrhopn;i++)
				fscanf(ff,"%lf %lf",&rl[i],&pl[i]);
			goto A;
			}
		fclose(ff);
		}
		
	ff=fopen("neutrino.dat","w");																			
	dlna=-log(amin)/(nrhopn-1);		
	printf("\n%e\n",amnu);
	fprintf(ff,"%e\n",amnu);	
													
	
	for(i=1;i<=nrhopn;i++)																					
	{																											
		a=amin*exp((i-1)*dlna);																				// an exponential scale has been taken for a						
		ninul(a,rhonu);																											
		rl[i]=log(rhonu[0]);																					//	store the density and pressure of the massive neutrinos as   							
		pl[i]=log(rhonu[1]);																					//	a function of a--(scale factor)											
		fprintf(ff,"%e  %e\n",rl[i],pl[i]);																													
		}			
	fclose(ff);																														
	
A:																																					
	splini();																									// this function should be called before calling splder										
	splder(rl,drl,nrhopn);																																		
	splder(pl,dpl,nrhopn);																																
	splder(drl,ddrl,nrhopn);																				// this will calculate the derivative of the array											
	
	fclose(ff);
	dq=1.0;
	for(i=1;i<=nqmax0;i++)
	{																											
		q=i-0.5;																									// this array is storing the values q^3 dq /(1+exp(q))	
		qdn[i]=dq*q*q*q/(exp(q)+1.0);																		// It will be used for calculating the neutrino energy	
		}																												
	}																												


/************************************************************************************************************************************/
/* This part is based on M.Zaldarriaga's Ph.D. Thesis																											*/
/* It will calculate the source term for the CMBR 																												*/
/************************************************************************************************************************************/

void foutput(int n,double y[],double yprime[],int j,double tau0,double tau,double zsanpass[])
{
	double x,a,a2,eta,etadot,alpha,alphadot,deltab,thetab;
	double thetabdot,deltag,thetagdot,sheargdot,polter,coupldot;
	double polterddot,shearrdot,rhonudot,shearnudot,omegavdyn;
	double weos,grhodot,grho,sgrhooa2,rgrho,adotoadot; 
	double dgsheardot,alphaddot,s1,s2,chi,chir,glens;
	double polterdot,ysanpass[2];  
	double d1,dp1,dk1,phi1;    

	d1=0.0;
	dp1=0.0;
	dk1=0.0;
	
	x=ak*(tau0-tau);
   a=y[1];
   a2=a*a;

   eta=y[3];
   etadot=yprime[3];

   alpha=(hdot+6*etadot)/(2.0*ak2);																				// M.Zaldarriaga's Ph.D. Thesis P.63
   alphadot=-3*dgshear/(2.0*ak2)+eta-2.0*adotoa*alpha;													// M.Zaldarriaga's Ph.D. Thesis P.63
   																														// Ma Bertschinger Eq(21d)

																															//  Baryons.
   deltab=y[6];
   thetab=y[7];
   thetabdot=yprime[7];
																															//  Photons.

   deltag=y[8];
   thetagdot=yprime[9];
   sheargdot=yprime[10]/2.0;

																															// Polarization term.
   polter=y[10]+y[9+lmaxg]+y[11+lmaxg];
   
   coupldot=8.0*(thetagdot+ak2*alphadot)/15.0;
   coupldot=coupldot-ak*0.6*(yprime[11]+yprime[10+lmaxg]+yprime[12+lmaxg]);
   
   polterdot=yprime[10]+yprime[9+lmaxg]+yprime[11+lmaxg];
   polterddot=coupldot-0.3*(dopac[j]*polter+opac[j]*polterdot);
  
																															// Massless neutrinos.
   shearrdot=yprime[12+2*lmaxg]/2.0;

																															// Second derivative of expansion rate
   if (amnu==0.0)																										// No massive neutrino
   {
	   rhonudot=0.0;
      shearnudot=0.0;
      }
   else																													// Massive neutrino present
   {
      nuder(a,rhonu,ysanpass,adotoa,y+iq2-1,yprime+iq2-1);
      rhonudot = ysanpass[0];
      shearnudot = ysanpass[1];
      }
      
																															// Dark energy
	omegavdyn = omegav*dynrho(a);
   weos = wdyn_func(a);

																															// Rho dot
   grhodot=(-grhom*(omegac+omegab)/a-2.0*(grhog+grhor*annur+grhonr*annunr*rhonu)/a2)*adotoa;
	grhodot=grhodot+grhonr*annunr*rhonudot/a2-((1.0+3.0*weos)*grhom*omegavdyn*a2)*adotoa;
	
	if(dimflag==1)
	{
		grho=pow((sqrt(grhom*(omegac+omegab)/a+(grhog+grhor*(annur+annunr*rhonu))/a2+grhom*omegav*a2)+sqrt(omegav*grhom)*a),2);
      sgrhooa2=sqrt(grho/a2);
      rgrho=sqrt(grhom*omegav+grhom*(omegab+omegac)/(a2*a)+(grhog+grhor*(annur+annunr*rhonu))/(a2*a2));
      grhodot=(sgrhooa2/rgrho)*(grhom*(omegab+omegac)*(-3/a)+(grhog+grhor*(annur+annunr*rhonu))*(-4/a2)+grhor*annunr*rhonudot/(a2*adotoa))+2*grho;
		}

	adotoadot=grhodot/(6*adotoa);																					// (ad/a)(d/dt)(ad/a) = (1/2)(d/dt)(ad/a)^2 
 	
																											//  Derivative of the shear
	dgsheardot=(4.0/3.0)*(grhog*sheargdot+annur*grhor*shearrdot)/a2;
	dgsheardot=dgsheardot-2.0*adotoa*dgshear+annunr*grhonr*shearnudot/a2;
	
																															// Calculation of the sources
	alphaddot=-3*dgsheardot/(2.0*ak2)+etadot-2.0*adotoadot*alpha-2.0*adotoa*alphadot;

   s1=etadot+alphaddot;																								// M.Zaldarriaga's Ph.D. Thesis Eq(3.38)
   s2=2*alphadot;																										// M.Zaldarriaga's Ph.D. Thesis Eq(3.38)

   d1=expmmu[j]*s1+vis[j]*(0.25*deltag+s2+polter/16+thetabdot/ak2+(3.0/16.0)*polterddot/ak2);
	d1=d1+dvis[j]*(alpha+thetab/ak2+3.0/8.0*polterdot/ak2)+ddvis[j]*(3.0/16.0)*polter/ak2;
	

   if (x>0.0)
   	dp1=vis[j]*(3.0/16.0)*polter/(x*x);																		//M.Zaldarriaga's Ph.D. Thesis Eq(3.38) Eq(3.37)
   else
      dp1=0.0;
        
																															// lensing visibility function; approximate epoch of recombination
   phi1=eta-adotoa*alpha;
   chi=tau0-tau;
   chir=tau0-taurmax;
   if (chi<chir) 																										// lensing convergence, expmmu is an approximation one should integrate 
   {																														// visibility function over r(chi)/r(tau) but the error is harmless
	   glens=(chir-chi)*chi/chir;
      dk1=glens*ak2*phi1*expmmu[j];
      }
   else
   	dk1=0.0;
   
   zsanpass[0]=d1;
   zsanpass[1]=dp1;
   zsanpass[2]=dk1;
   zsanpass[3]=phi1;	
	}

/************************************************************************************************************************************/  
/*		It will store the transfer function at different point																								*/
/*		It will be used by the function cmbfalt()																													*/
/************************************************************************************************************************************/  

void outtransf(int n,double y[],double curv,int itf)
{
	double tfc,tfb,tfg,tfn,a,deltan,tfnm,beta,pnu,rhonu,xsanpass[2];
	double zsanpass[4];
	
	beta=ak;//sqrt(ak-curv);
   ak2=beta*beta;

   tfc=y[4]/ak2;
   tfb=y[6]/ak2;
   tfg=y[8]/ak2;
   tfn=y[10+2*lmaxg]/ak2;

   if(amnu!=0.0)
   {
	   a=y[1];
      nu1(a,xsanpass);
      rhonu=xsanpass[0];
      pnu=xsanpass[1];
      nu2(a,zsanpass,y+iq0,y+iq1,y+iq2);
      drhonu = zsanpass[0];
      
      deltan=drhonu/rhonu;
      tfnm=deltan/ak2;
      }

	power[6*istore[itf]+1][itf] =beta/(h0/100.0);
   power[6*istore[itf]+2][itf] =tfc;
   power[6*istore[itf]+3][itf] =tfb;
   power[6*istore[itf]+4][itf] =tfg;
   power[6*istore[itf]+5][itf] =tfn;

   if (amnu==0.0)
   	power[6*istore[itf]+6][itf] =0.0;
   else
   	power[6*istore[itf]+6][itf] =tfnm;
      
   istore[itf] = istore[itf] +1;
   
	}


/************************************************************************************************************************************/
/*		This function will initiallize the tensor perturbation terms																						*/
/************************************************************************************************************************************/

void finitialt(double tau,double y[])
{
	int ind1,ind2,ind3;
	double a;
		
   a=tau*adotrad;
   y[1]=a;

																											// Tensor modes
   y[2]=1.0;																							// ht=1.0
   y[3]=0.0;																							// htpr=0.0
   ind1=4;
   ind2=ind1+lmaxt+1;
   ind3=ind2+lmaxt+1;
   for(int l=0;l<=lmaxt;l++)
   {
   	y[ind1+l]=0.0;
   	y[ind2+l]=0.0;
   	}
   for(int l=0;l<=lmaxnr0;l++)
   	y[ind3+l]=0.0;
   }


/************************************************************************************************************************************/     
/*  This function will calculate the time derivatives of the tensor perturbation																		*/
/*  This part is based on the followingg papers :																												*/
/*																																											 	*/	
/*	 1) A Beginner’s Guide to the Theory of CMB Temperature and Polarization Power Spectra in the Line-of-Sight Formalism            */
/*		 :: Yen-Ting Lin and Benjamin D. Wandelt  																												*/	
/*		 :: Represented by : Y                                                              														*/
/*	 2) CMB Anisotropies: Total Angular Momentum Method :: Wayne Hu & Martin White																	*/
/*	 	 :: Represented by : W																																			*/	
/*	 3) Cosmic Microwave Background Polarization :: Kosowsky																									*/
/*	 	 :: Represented by : K																																			*/
/*	 4) Imprint of Gravitational Waves on the Cosmic Microwave Background 																				*/
/*	    :: Robert Crittenden, J.Richard Bond, Richard L.Davis, Gorge Efstathiou, Paul J. Steinhardt												*/
/*	    :: Represented by : R	 	  																																	*/		
/*																																											 	*/			
/************************************************************************************************************************************/  

void fderivst(int n,double x,double y[],double yprime[])
{
	double ep,tau,a,a2,tcp,tcp1,tcp2,rhonu,pnu,omegavdyn,grho;
	double shearg,shearr,pi,ht,htpr,psie,deltat0,deltap0;
	double ep0=1.0e-2,xsanpass[3],htdpr,opac1,ysanpass[2];
	int ind1,ind2,ind3;

   if (ak>0.06*epsw)																					// ep is used to stop the tight coupling approximation.
   	ep=ep0;
   else
   	ep=1.17*ep0;																													

	tau=x;
   a=y[1];																								// Time of the calculation		

   a2=a*a;
   
   thermo(tau,xsanpass);																			// Calculate Opacity
   opac1=xsanpass[1];

																											// Tight Coupling parameters
   tcp=0.0;
   tcp1=ak/opac1;
   tcp2=1.0/(opac1*tau);
   if ((tcp1>ep)||(tcp2>ep))
   	tcp=1.0;																							// Not tightly Coupled 
      
																											// Compute expansion rate.
   if (amnu==0)																						// No massive Neutronos
   {
   	rhonu=1.0;
      pnu=1.0/3.0;
      }		
   else
   	nu1(a,ysanpass);																				// Masssive Neutrinos Present, So calculate pressure and density
   rhonu = ysanpass[0];
   pnu = ysanpass[1];
        
   omegavdyn = omegav*dynrho(a);
   
																											//  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
   grho=grhom*(omegac+omegab)/a+(grhog+grhor*annur+grhonr*annunr*rhonu)/a2+grhom*omegavdyn*a2;

	if(dimflag==1)																						// Energy density for  5 dimensions
	{
		grho=pow((sqrt(grho)+sqrt(omegav*grhom*a2)),2);
		}

   adotoa=sqrt(grho/3.0);
   yprime[1]=adotoa*a;																				 	

   ind1=4;
   ind2=ind1+lmaxt+1;
   ind3=ind2+lmaxt+1;
   shearg=y[ind1]/15+y[ind1+2]/21+y[ind1+4]/35;
   shearr=y[ind3]/15+y[ind3+2]/21+y[ind3+4]/35;

	/*********************************************************************************************************************************/  
	/*  Relativistic components																																		*/
	/*********************************************************************************************************************************/  
   pi=(grhog*shearg+annur*grhor*shearr)/(a2*3);												// 

	/*********************************************************************************************************************************/  
	/*  Tensors																																								*/
	/*********************************************************************************************************************************/  

   ht=y[2];
   htpr=y[3];
   yprime[2]=htpr;
   htdpr=-2*adotoa*htpr-ak2*ht+24*pi;															//Eq(W.73), Eq(R.3)
   yprime[3]=htdpr;																						

	/*********************************************************************************************************************************/  
	/*  Photon perturbations																																			*/			
	/*********************************************************************************************************************************/  
	psie=y[ind1]/10.0+y[ind1+2]/7.0+3*y[ind1+4]/70.0-3.0*y[ind2]/5.0+6.0*y[ind2+2]/7.0-3.0*y[ind2+4]/70.0;
																											//Eq(Y.61), Eq(R.3)
//	printf("tau = %e psie %e opac %e\n",tau,psie,opac1);																														
   if (tcp==1) 																						// Not Tightly Couppled	
   {
	   yprime[ind1]=-ak*y[ind1+1]-opac1*y[ind1]+opac1*psie-htpr;							// Eq(R.3)
      yprime[ind2]=-ak*y[ind2+1]-opac1*y[ind2]-opac1*psie;									// Eq(R.3)

      for(int l=1;l<=lmaxt-1;l++)
      {
	      yprime[ind1+l]=ak*(l*y[ind1-1+l]-(l+1)*y[ind1+1+l])/(2*l+1)-opac1*y[ind1+l]; // Ma, Bertschinger Eq(64)
         yprime[ind2+l]=ak*(l*y[ind2-1+l]-(l+1)*y[ind2+1+l])/(2*l+1)-opac1*y[ind2+l]; // Ma, Bertschinger Eq(64)
         }

																											// Truncate moment expansion
																											// Ma, Bertschinger Eq(64)
      yprime[ind1+lmaxt]=ak*y[ind1-1+lmaxt]-(lmaxt+1)/tau*y[ind1+lmaxt]-opac1*y[ind1+lmaxt];
      yprime[ind2+lmaxt]=ak*y[ind2-1+lmaxt]-(lmaxt+1)/tau*y[ind2+lmaxt]-opac1*y[ind2+lmaxt];
      }  

	else																									// Tightly Couppled
   {
	   deltat0=-4.0*htpr/(opac1*3.0);
      deltap0=-deltat0/4.0;
      y[ind1]=deltat0;
      y[ind2]=deltap0;
      
      for(int l=0;l<=lmaxt;l++)
      {
	      yprime[ind1+l]=0.0;
         yprime[ind2+l]=0.0;
         }
 		}
        
	/*********************************************************************************************************************************/  
	/*  Massless Neutrino perturbations																																*/			
	/*********************************************************************************************************************************/  
	yprime[ind3]=-ak*y[ind3+1]-htpr;

   for(int l=1;l<=lmaxnr-1;l++)
   	yprime[ind3+l]=ak*(l*y[ind3-1+l]-(l+1)*y[ind3+1+l])/(2*l+1);					// Ma, Bertschinger Eq(49)

																											// Truncate moment expansion
	yprime[ind3+lmaxnr]=ak*y[ind3-1+lmaxnr]-(lmaxnr+1)/tau*y[ind3+lmaxnr];			// Ma, Bertschinger Eq(49), Eq(51) Combined
//	for(int ab=1;ab<=n;ab++)
//		printf("%e\t%e\n",y[ab],yprime[ab]);
//	printf("\n");
//	exit(1);	
	}


/************************************************************************************************************************************/
/* This part is based on 																																				*/
/*	1) Signatures of a Graviton Mass in the Cosmic Microwave Background																					*/
/* 	:: Sergei Dubovsky, Raphael Flauger, Alexei Starobinsky, Igor Tkachev																			*/
/*	2)	Microwave Background Constraints on Cosmological Parameters																							*/
/*		:: Matias Zaldarriaga, David N. Spergel, Uro¡s Seljak																									*/
/*	3) Signature of Gravity Waves in Polarization of the Microwave Background 																			*/
/*		:: Uro¡s Seljak, Matias Zaldarriaga																															*/
/* It will calculate the tensor source term for the CMBR 																									*/
/************************************************************************************************************************************/

void foutputt(int n,double y[],double ypr[],int j,double tau0,double tau,double xsanpass[])
{
	double x,x2,htpr,htdpr,psie,psiedot,psieddot;
	double dt1,dte1,dtb1;
	int ind1,ind2;
	
	x=ak*(tau0-tau);
   x2=x*x;
   ind1=4;
   ind2=ind1+lmaxt+1;
   htpr=y[3];
   htdpr=ypr[3];
   
  																															// Sergei Eq(8), Zaldarriaga Eq(6)
   psie=y[ind1]/10.0+y[ind1+2]/7.0+3.0*y[ind1+4]/70.0-3.0*y[ind2]/5.0+6.0*y[ind2+2]/7.0-3.0*y[ind2+4]/70.0;
   psiedot=ypr[ind1]/10.0+ypr[ind1+2]/7.0+3.0*ypr[ind1+4]/70.0-3.0*ypr[ind2]/5.0+6.0*ypr[ind2+2]/7.0-3.0*ypr[ind2+4]/70.0;
   psieddot=-0.3*(opac[j]*psiedot+dopac[j]*psie)-0.1*htdpr-ak*(3.0*ypr[ind1+1]/70.0+ypr[ind1+3]/15.0+ypr[ind1+5]/42.0-33.0*ypr[ind2+1]/35.0+8.0*ypr[ind2+3]/15.0-ypr[ind2+5]/42.0);

   if (x>0.0)
   {
	   dt1=(-expmmu[j]*htpr+vis[j]*psie)/x2;																	// Zaldarriaga Eq(6)
      dte1=vis[j]*(psie-psieddot/ak2-6.0*psie/x2-4.0*psiedot/ak/x)-dvis[j]*(4.0*psie/x/ak+2.0*psiedot/ak2)-ddvis[j]*psie/ak2;
      dtb1=2.0*(vis[j]*(2.0*psie/x+psiedot/ak)+dvis[j]*psie/ak);										// Zaldarriaga Eq(6)
      }
   else
   {
	   dt1=0.0;
      dte1=0.0;
      dtb1=0.0;
      }
   
   dte1=-dte1;
   dtb1=-dtb1;
   
   xsanpass[0] = dt1;
   xsanpass[1] = dte1;
   xsanpass[2] = dtb1;
   }


/*************************************************************************************************************/
/* 	This function will short the k arrays in the increasing order														 */
/*************************************************************************************************************/

int indexx1(int n,double arr[],int indx[])
{
	int jstack,l,ir,indxt,i,itemp,k;
	double asan;
	int M=7,NSTACK=50;
	int istack[NSTACK+1];
	int j;
//	int i;
	for(int js=1;j<=n;j++)
   	indx[js]=js;

	jstack=0;
   l=1;
   ir=n;
//   printf("%d\n",n);
//   for(int ii=1;ii<=n;ii++)
//  		printf("%e %d\n",arr[ii],indx[ii]);
 //  exit(1);
   
   while(1)
   {
	if(ir-l<M)
	{
//		printf("ir,l = %d %d",ir,l);
		for(int j=l+1;j<=ir;j++)
		{
			indxt=indx[j];
         asan=arr[indxt];
         i=j-1;
         for(i=j-1;i>=1;i--)
         {
	         if(arr[indx[i]] <= asan) break;
            indx[i+1]=indx[i];
//            printf("%d %d\n",i,j-1);
				}
//         printf("Hi %d",i);
         if(arr[indx[i]] > asan)
         	i=0;
//         printf("Cought"); 	
         indx[i+1]=indxt;
			}
			
//		printf("ir,l = %d %d",ir,l);
				        
      if(jstack==0)
      	return 0;
      
		ir=istack[jstack];
      l=istack[jstack-1];
      jstack=jstack-2;
// 		printf("%d %d %d\n",ir,l,jstack);
      }
      
        
   else
   {
      k=(l+ir)/2;
      itemp=indx[k];
      indx[k]=indx[l+1];
      indx[l+1]=itemp;

		if(arr[indx[l+1]] > arr[indx[ir]])
      {
	      itemp=indx[l+1];
         indx[l+1]=indx[ir];
         indx[ir]=itemp;
        	}
        		
      if(arr[indx[l]] > arr[indx[ir]])
      {
	      itemp=indx[l];
         indx[l]=indx[ir];
         indx[ir]=itemp;
        	}
      
      if(arr[indx[l+1]] > arr[indx[l]])
      {
	      itemp=indx[l+1];
         indx[l+1]=indx[l];
         indx[l]=itemp;
        	}
      
      i=l+1;
      j=ir;
      indxt=indx[l];
      asan=arr[indxt];

SAN3:
      i=i+1;
      
      if(arr[indx[i]] < asan) goto SAN3;

		do
		{
			j=j-1;
			}while(arr[indx[j]] > asan);
      
      if(j >= i)
      {
      
	   	itemp=indx[i];
      	indx[i]=indx[j];
      	indx[j]=itemp;
     		goto SAN3;
      	}
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack=jstack+2;
      
      if(jstack > NSTACK)
      {
	      printf("NSTACK too small in indexx");
	      exit(1);
	      }
      
      if(ir-i+1 >= j-l)
      {
			istack[jstack]=ir;
         istack[jstack-1]=i;
         ir=j-1;
         }
      else
      {
	      istack[jstack]=j-1;
         istack[jstack-1]=l;
         l=i;
         }
   	} 
   	}  
   	return 0;	
	}


/******************************************************************************************************************/
/*		It will printout the transfre function at the given points.																	*/
/******************************************************************************************************************/

void output_power(int ntf,double amnu)
{	
	int kmax1=1000;
	int nmax=1000;
	int indx[1001]={0};
		
	double k[kmax1+1],tfc[kmax1+1],tfb[kmax1+1],tfg[kmax1+1],tfn[kmax1+1],tfnm[kmax1+1];
	FILE *fp;
	fp = fopen("transfcmb.d","w");
	for(int itf=1;itf<=ntf;itf++)
	{
   	if(istore[itf] > nmax)
   	{
      	printf("\n Need to increase nmax");
         exit(1);
         }
      for(int ncount=1;ncount<=istore[itf];ncount++)
      {
      	k[ncount] = power[6*(ncount-1)+1][itf];
         tfc[ncount] = power[6*(ncount-1)+2][itf];
         tfb[ncount] = power[6*(ncount-1)+3][itf];
         tfg[ncount] = power[6*(ncount-1)+4][itf];
         tfn[ncount] = power[6*(ncount-1)+5][itf];
         tfnm[ncount] = power[6*(ncount-1)+6][itf];
         }

		/*************************************************************************************************************/
		/*  	sort k's in increasing order																									 */
		/*************************************************************************************************************/

	   indexx(istore[itf],k,indx);
      for(int icount=1;icount<=istore[itf];icount++)
      {
      	if(k[indx[icount]]!=0.0)
      	{
         	if (amnu!=0.0)
            	fprintf(fp,"%e %e %e %e %e %e\n",k[indx[icount]],tfc[indx[icount]],tfb[indx[icount]],tfg[indx[icount]],tfn[indx[icount]],tfnm[indx[icount]]);
            else
            	fprintf(fp,"%e %e %e %e %e\n",k[indx[icount]],tfc[indx[icount]],tfb[indx[icount]],tfg[indx[icount]],tfn[indx[icount]]);
            }
         }
		}
	fclose(fp);	
	}


/*************************************************************************************************************/
/*     This subroutine reads the jl files from disk and initializes other variables needed in CMBANUS.		 */
/*************************************************************************************************************/

void initjl()
{
	int kmaxfile,lmofile,l0file;//,kmax;
	int ketamax=3*l0max+126,lfile[lmax+1];
	double xlimmin=10.0;
	double x,xlim;
	double d0hi=1.0e40,d0lo=1.0e40;
	double aj=0.0;
	FILE *fp,*fp34;
	
	fp = fopen("jlgen.dat","r");
   fp34=fopen("newdat.d","w");
   fscanf(fp,"%d %d",&lmofile,&kmaxfile);
																									// The code needs 300 more ls for the lensing calculation.
	if ((lmo > lmofile)||(akmax0 > kmaxfile))
   {
	   printf("\n\nThe value of lmax and/or kmax is more than what is stored in the file");
	   printf("\nIn the file we have lmax=%d and kmax = %d",lmofile-300,kmaxfile);
      printf("\nSo restart the calculation please. \n");
      exit(1);
      }

	kmaxfile=kmaxfile-25+151;
   kmaxjl=(akmax0-25+151);

	double bjl[kmaxfile+1],bjlpr[kmaxfile+1];   

   if (kmaxfile > ketamax)
   {
	   printf("kmax in file %d is too large. Maximum value for ketamax can be %d", kmaxfile,ketamax);
      exit(1);
      }


	fscanf(fp,"%d",&l0file);																// Checking if the lvalues.inc file used to build jl file is the same as the one in the code.
   for(int j=1;j<=l0;j++)
   {
   	fscanf(fp,"%d",&lfile[j]);
   	printf("%d\n",lfile[j]);
   	}

	for(int j=1;j<=l0;j++)
   {
	   if (l[j]!=lfile[j])
	   {
		   printf("\n\nlvalues.inc file used to build jl file and the one in the code differ.");
         printf("\nYou must use the same one %d %d %d",j,l[j],lfile[j]);
         exit(1);
         }	
      }
																								
	printf("\nkmaxfile = %d\n",kmaxfile);
   for(int i=1;i<=kmaxfile;i++)														// reading  j_l remember to create jl.dat with jlgen.f first using correct lmax and ketamax
   {
		if (i<151)
	   {
		   if (i<=51)
         	xx[i]=double(i-1)/10.0;
         else
         	xx[i]=double(i-51)/5.0+5.0;
         }
      else
      	xx[i]=(i-151)+25.0;
      }
   printf("Hi %d",l0);   
   for(int j=1;j<=l0file;j++)
	{
		for(int i=1;i<=kmaxfile;i++)
	   {
		   x=xx[i] ;
         xlim=(l[j])/20.0;
         xlim=max(xlim,xlimmin);
			
         xlim=l[j]-xlim;
         if (x-xlim > 1e-10)
         {
	         aj=0.0;
	         fscanf(fp,"%lf",&aj);
            ajl[i][j]=aj;
            }
         else
         	ajl[i][j]=0.0;
 			}
 		}
   fclose(fp);
	fclose(fp34);
   for(int i=2;i<=(kmaxfile-1);i++)
   {
	   dxx[i]=(xx[i+1]-xx[i-1])/2.0;
 		}
	  		
   dxx[1]=xx[2]/2.0;
   dxx[kmaxfile]=(xx[kmaxfile]-xx[kmaxfile-1])/2.0;

   bjl[0]=0;
   bjlpr[0]=0;


   for(int j=1;j<=l0;j++)																	// Get the interpolation matrix for bessel functions
   {
	   for(int k=1;k<=kmaxfile;k++)
	   {
		   bjl[k]=ajl[k][j];
 			}	   	
	   spline(xx,bjl,kmaxfile,d0lo,d0hi,bjlpr);
	   
	   for(int k=1;k<=kmaxfile;k++)
	   	ajlpr[k][j]=bjlpr[k];
//	   	printf("%e\n",bjlpr[k]);
 		}
//	exit(1);
   }



 
	

/*************************************************************************************************************/
/* 	This subroutine computes the power spectra for mode ak of the tensor perturbations.							 */	
/*************************************************************************************************************/

void powertflat(double ak,int in,double *apower)
{
	double aknlog;
	aknlog=log(ak/0.05);
   *apower=exp((ant[in]+.5*alphant[in]*aknlog)*aknlog);
	}


/************************************************************************************************************************************/
/*		All the calculation of the CMBFLAT is done in this function. This is the main function of this program. Some of the Cl will  	*/
/*		be calculated in  this program other Cl's will be manupulated by fitting an spline															*/
/************************************************************************************************************************************/

void testveriables(int nvar,double tau,double y[],double tauend,double yprime[])
{
	printf("tau = %e tauend = %e\n",tau,tauend);
	for(int i=1;i<=nvar;i++)
	{
		printf("y[%d] = %.8e\typrime[%d] = %.8e\n",i,y[i],i,yprime[i]);
		}
	exit(1);	
	}

void findis(int indtau[])
{
	printf("ntfmax = %d ntau = %d",ntfmax,ntau);
	printf("\nIstore: \n");
	for(int i=1;i<=ntfmax+1;i++)
		printf("%d\t",istore[i]);
	printf("\nIndtau: \n");
	for(int i=1;i<=ntau+ntfmax+1;i++)	
		printf("%d\t",indtau[i]);
//	printf("\n%d",indtau[ntau+ntfmax+1]);	
	exit(1);
	}
	

void cmbflat()
{
	double omegam;																								// Matter density
	double zst;																									// Redshift upto which the equations will be integrated
	double arel=0.0;																								
	double akmax,akmin,dlnk;																				// akmin -> k from which the integration will start
																													// akmax -> integration will be done till this k value
																															
	double taumax,taumin,dlntau0,taurend;																// taumax -> where to stop intigration
																													// taumin -> where to start intigration	
																													// 

	double tautf[ntfmax+1];																					// tautf <- stores the conformal time at which the transfer functions are requested
	
	
	double shor;																								//sound horizon			
	
	int n1;
	int nstep;																									// Number of discritization of the time step

	double atf;																									// Tempurary variable	
//	double tauf[ntfmax+1];
	
	double tautf0[ntau+ntfmax+1];																			// Discritized the time																						
	int indtau[ntau+ntfmax+1]={0};																				// It will store where we have to calculate the transfer functions
	indtau[ntau+ntfmax+1]=0;
	
	double dkn1,dkn2,dlnk0;
	double tau,taunew=0;
	int nriend,nk,nk1,nk2,n10;
	int kkmin, kkmax;
	
	double akchange;
	double ak0[nk0+1],ak10,t10,a0,ar;																	// ak0[] <- Store the wavenumbers at the gridpoints	
	double tauend,taulasttf,taustart,taustart0;
	int itau,ind;																								// itau is a flag
	
	double phi1,akt,ho,xf,b0;//,bo
	int closestknum;
	double tol1;																								// tolerance vaue
	
	double y[nvar0+1],yprime[nvar0+1],yt[nvar0t+1],ytprime[nvar0t+1];
   double zsanpass[4],phi,stpt,ysanpass[3],akdone,akgf,phihk;
   	
   double c[25];																								// dvark intigrator communication parameter
	double gftemp1,gftemp2,gftemp3,gftemp4,gftemp5,gftemp6,growthfactor;
	double d0lo=1.0e40,d0hi=1.0e40;
	double dk,dk0,dlnk1;
	double ak1[nkmax+1],dak1[nkmax+1];
	
   double dtl[lmax+1],dtl2[lmax+1],dtl3[lmax+1];
   double dtel[lmax+1],dtel2[lmax+1],dtel3[lmax+1];
   double dtbl[lmax+1],dtbl2[lmax+1],dtbl3[lmax+1];
   double st2[nstep0+1],ste2[nstep0+1];
   double stb2[nstep0+1];

   double dl[lmax+1]={0},dl2[lmax+1]={0},dl3[lmax+1]={0};
   double dpl[lmax+1]={0},dpl2[lmax+1]={0},dpl3[lmax+1]={0};
   double dkl[lmax+1]={0},dkl2[lmax+1]={0},dkl3[lmax+1]={0};
   double s2[nstep0+1],sp2[nstep0+1],sk2[nstep0+1];	

   double dladb[lmax+1]={0},dl2adb[lmax+1]={0},dl3adb[lmax+1]={0};
   double dpladb[lmax+1]={0},dpl2adb[lmax+1]={0},dpl3adb[lmax+1]={0};
   double dkladb[lmax+1]={0},dkl2adb[lmax+1]={0},dkl3adb[lmax+1]={0};
   double s2adb[nstep0+1],sp2adb[nstep0+1],sk2adb[nstep0+1];	

   double dliso[lmax+1]={0},dl2iso[lmax+1]={0},dl3iso[lmax+1]={0};
   double dpliso[lmax+1]={0},dpl2iso[lmax+1]={0},dpl3iso[lmax+1]={0};
   double dkliso[lmax+1]={0},dkl2iso[lmax+1]={0},dkl3iso[lmax+1]={0};
   double s2iso[nstep0+1],sp2iso[nstep0+1],sk2iso[nstep0+1];	

   double de2,xlim,xlmax1,xlmax2,tmin,tmax,h2,a2,b2;     
   double ddt1,ddt2,ddp1,ddp2,ddk1,ddk2,ddt,ddp,ddk;
   double ddt1adb,ddt2adb,ddp1adb,ddp2adb,ddk1adb,ddk2adb,ddtadb,ddpadb,ddkadb;
   double ddt1iso,ddt2iso,ddp1iso,ddp2iso,ddk1iso,ddk2iso,ddtiso,ddpiso,ddkiso;

   double xlimmin=10.0;  
   double dtau3,xi,x;
   double ajl0;
   double ckj,cpkj,cckj,ckkkj,ctkkj;
   double dtdt1,dtdt2,dtde1,dtdb1,dtde2,dtdb2,dtdt,dtdb,dtde;
   double ctkj,ctekj,ctbkj,ctckj,apowert;
   double cllo,clhi;
   double ctnorm,xl[lmax];
   double clint,cplint,cclint,ctklint; //,ckklont
	double clintadb,cplintadb,cclintadb;
	double clintiso,cplintiso,cclintiso;
	double clintcross,cplintcross,cclintcross;
	double clinttotal,cplinttotal,cclinttotal;
   double ctlint,ctelint,ctblint,ctclint,ckklint;
   
	int nstart1,nstop1,nstop1a,nstart2,m2,nstop2,nstop2a;
	int no,no1,nko;
	int nkt,nstps,nstpt;
	int khi,klo;
	int mxx[nstep0+1],m1;
	int llo,lhi; 
	int nw,nwt;
	double w[10][nvar0+1],wt[10][nvar0+1];

	nw = nvar0;
	
   dl2[lmax]=0;
   dl3[lmax]=0;
   dpl2[lmax]=0;
   dpl3[lmax]=0;
   dkl2[lmax]=0;
   dkl3[lmax]=0;


	omegam = omegab+omegan+omegac;																		// Total matter density of the universe Baryonic + Neutrino + CDM
	
	// (P1) : Calculating the final redshift at which the program will stop the calculations 
																													// arXiv:astro-ph/9603033					
																													// Seljak, Zaldarriaga (1996)  -  Page No 13. "Free streaming"
																																
	if(fabs(omegab+omegac-1.0)<.001 && zri==0.0 && optdlss==0.0 && itflag==0 && h0>40.0) // The poin where the program stop calculation. If reionization or 
		zst=10.0;																								//	the tensor perturbations are requested then we should intigrate
	else 																											// till the present era	
		zst=0.0;																	
		
	if(ict!=0)																									// If transfer functions are requested then we should stop at the z 	
		zst=min(zst,ztf[ntf]);																				// where the last transfer function is requested or after that 
	
	// (P2) : 
		
	if(itflag!=2) 																								// If scalar terms are not present then photons equations will be  
		lmaxg=lmax0;																							// discritized into these many steps
	
	lmaxnr=lm2;																							      //	Number of discritization of the relativistic Neutrino	
	nstep=nstep0;																								//	Number of discritization of the time step
	

	if(itflag !=0)																								// If tensor perturbations are requested	
		lmaxt=lmaxt0;																							// Number of discritization of the tensor equations.
	else
		lmaxt=0;

	/*********************************************************************************************************************************/
	/*		Calculate number of equations																																*/
	/*																																											*/
	/*********************************************************************************************************************************/
    
   if (itflag!=2)																								// Scalers are requested
   {
	   iq0=11+2*lmaxg+lmaxnr;																				// iq0, iq1, iq2 will be used as an indicator  
      iq1=iq0+nqmax;																							// when we will store the massive neutrino number density,
      iq2=iq1+nqmax;																							//	density and pressure, in function nu2()
      
      if (ndyn==0 || ndyn==3 || ndyn==4)																/********************************************************************/
      { 																											/*  1) 3 				for storing a, adot, eta							  */
	      nvar=7+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1);											/*	 2) 2 				for storing CDM										  */
         }																										/*	 3) 2 				for storing Baryon									  */
      																											/*  4) 2 				for photon intensity and polarisation			  */
      else																										/*  5) 2*lmaxg 		for storing photons shears							  */
      {																											/*	 6) lmaxnr  		for storing massless neutrino 					  */
	      nvar=9+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1);											/*	 7) 3*nqmax			for storing massive neutrino						  */
         }																										/*	 8) nqmax*lmaxnu 	for Psi_nu for all momentum of massive neu	  */
      }																											/*  9) 2 				for dark energy perturbation						  */	
      																											/*	 																			  		  */
      																											/*	  This is how the above numbers came :)								  */																		  
   else 																											/********************************************************************/								
   {
	   nvar=0;																									// nvar is for dark energy. It will be used in function finitial(),fderivs() 
      }
      
   if (itflag!=0)																								// Tensors are requested
   {
	   nvart=(2*lmaxt+5)+lmaxnr+1;																		/*********************************************************************/																
	   }																											/*	 1) 3					for a, h, h_dot											*/
	else																											/*	 2) 2*lmaxt			for storing h,h_dot									 	*/
	{																												/*  3) lmaxnr			for massless neutrinos									*/		
		nvart=0;																									/*																							*/
		}																										   /*		From this we get the following numbers	:)								*/
													                                                /*********************************************************************/                     																											

   // (P3) : Massive neutrino distribution 

	/****************************************************************************************************************/
	/*		Initiallize Massive Neutrino distribution										  										    	 */
	/*		"Cosmology"	- S.Weinberg																											 */
	/****************************************************************************************************************/
	
	if(annunr==0 && omegan==0.0)																			// If no massive neutrinos are present then set  
	{																												// all the massive neutrino parameters to zero	
		amnu=0.0;																											
		nqmax=0;																															
	   lmaxnu=0;																												
		}																																		
																															
	else																											// If neutrino present then
	{																															
	   amnu=(omegan/annunr)*(grhom*1000/grhonr)*constan/(1.5*zeta3);							// amnu=m_nu*c**2/(k_B*T_nu0)
		nqmax=nql;
		lmaxnu=lm3;
		
		double q,dq;																							// Tempurary variables :)
		dq=1.0;

		for(int i=1;i<=nqmax;i++)
		{
			q=i*dq-0.50;
			dlfdlq[i]=-q/(1+exp(-q));																		// Calculate the d(log f)/d(log q)
			}
		}
	
	/****************************************************************************************************************/
	/*		Initiallize Massive Neutrino distribution																				 		 */
	/****************************************************************************************************************/
	
	if((annunr!=0)||(omegan!=0.0))
	{
		arel=1.0-3.0/amnu;
		initnul();
		}
	
		
	// (P4) : Calculating time of different epochs	
	
	/****************************************************************************************************************/
	/*		Time at the present epoch																											 */
	/****************************************************************************************************************/	

	tau0=rombint(dtauda,1.0e-8,1.0,tol);															// This gives the conformal time in the present era
	epsw=100.0/tau0;																						// This will set a limit on the wavenumber during the intigration process

	/*****************************************************************************************************************/	
	/*		Calculate the time of reionization																								  */
	/*****************************************************************************************************************/	
	
	// (P5) : Calculating reionization
	
	if(optdlss>0.0)
      reiopar();																							// It will calculate zri and zristp
   
   // (P6) :  k and t grid speacification
   
   akmax=akmax0/tau0;																					// Maximum and minimum k-values. 
   akmin=0.15/tau0;
   
   if (itflag==0)
   	dlnk=0.1;																							// dlnk is the logarithmic spacing used for low k.
   else																										// if only scalars are requested then we large speacing can be taken		
   	dlnk=0.05;
  
   if(akmax==0)
   	dtaurec=0.0;																						// This will only happen if there is some error in the user input for "ketamax" 
   else
   	dtaurec=4.0/akmax;   
   
   taumax=rombint(dtauda,1.0e-8,1.0/(zst+1.0),tol);											// taumax is the time where the program stops calculation

   dlntau0=0.0050;

   if (itflag!=0) dlntau0=0.0025;																	// if only tensor perturbations are required then take small step 	
     
   if (zri!=0.0)
   {
      taurist=rombint(dtauda,1.0e-8,1.0/(1+zri),tol);											// Conformal time where the reionozation starts 		
      tauristp=rombint(dtauda,1.0e-8,1.0/(1+zristp),tol);									// Conformal tile where the reionization stops
      }
   else
   {
	   taurist=tau0;
      tauristp=tau0;
      }
   	
   if(ict==0) taumin=0.001/akmax;																	// taumin if only Cl is requested 
   if(ict==1) taumin=0.001/akmaxt;																	// taumin if only TF is requested
   if(ict==2) taumin=0.001*min(1.0/akmax,1.0/akmaxt);											//	taumin if both are needed 
	
	taumin=min(taumin,0.1);																				// Start in the radiation dominated era 
	  
   if(amnu!=0.0)
   	taumin = min(taumin,arel/adotrad);															// Start when massive neutrinos were strongly relativistic
   
	/*************************************************************************************************************/
	/* 			Initialize baryon temperature and ionization fractions vs. time.											 */
	/*************************************************************************************************************/   
   double tausan[2],dlsso;																				
   int nstepabc[1];																						
   nstepabc[0] = nstep;																					// 
	n1 = finithermo(taumin,taumax,dlntau0,nstepabc,tausan);									//output - > taurend, n1,
	nstep = nstepabc[0];																														
	taurend = tausan[0];																								
	
	// (P7) : k splitting 
																														
	dlsso = tau0-taurmax;																							
   shor=rombint(dsoundda,1.0e-8,armax,tol);																// Sound Horizon at LSS

	aksplit=(1.5/shor)*aksplit;																				// K splitting
																														// Compute cut wavevector
  																														
   if ((kcutflag==1)||(kcutflag==-1))																		// Output dlss to shift spectra.
   {
	   printf("### First line has dlss for this model and kcutflag");								// Write dlss in the output file.
	   printf("### Then normal l ClT ClE ClC output");
      printf("%e %d",dlsso,kcutflag);
      }
   
	// (P8) : Calculating the points where the transfer functions will be stored
      
	int itf;
	double aprevious;
	if (ict!=0)																										// Transfer functions are required	
	{
		a0=1.0e-8;
		tautf[0] = 0.0;
		aprevious = a0;
		
	   for(itf=1;itf<=ntf;itf++)																				// Calculating the times for the outputs of the transfer functions.
	   {
  	 		atf=1.0/(ztf[itf]+1.0);
      	tautf[itf]=tautf[itf-1]+rombint(dtauda,aprevious,atf,tol);
      	tautf[itf]=min(tautf[itf],atau0[nstep]);
      	indtau[itf]=itf;
      	aprevious = atf;
			}	
		
		if(lensflag!=0)
   	{
		   for(itf=1;itf<=ntf;itf++)
      		tautf0[itf]=tautf[itf];
      	
      	for(int i=1;i<=ntau+ntf;i++)
	    		indtau[i]=0;
      	
      	ar=1.0/1090.0;
      	taur=rombint(dtauda,a0,ar,tol);
      	itf=1;
      	
      	taunew=taur;
      	
      	for(int i=1;i<=ntau;i++)
      	{
         	if (taunew<=tautf0[itf] || itf>ntf)
         		tautf[i+itf-1]=taunew;
         
         	else
         	{
	         	tautf[i+itf-1]=tautf0[itf];
            	indtau[i+itf-1]=itf;
            	itf=itf+1;
            	continue;
            	}
         	taunew=i*(tau0-taur)/(ntau-1)+taur;
         	}
						            
			if(itf<=ntf)
      	{
	      	for(int i=itf;i<=ntf;itf++)
         	{
	         	indtau[ntau+i]=i;
            	tautf[ntau+i]=tautf0[i];
            	}
         	}
      	ntf=ntf+ntau;
      	}
   	}
 
		
   tautf[ntf-1]=tautf[ntf-1]-tol;       
																						// Integration will only be carried out after z=10 for low
																						// k, k < k10. If there is reionization the boundary will not be z=10 but tauristp
	// (P9): No need to calculate for the large wavelengths after z=10

   if (zst!=10.0)
   {
   	t10=rombint(dtauda,1.0e-8,1.0/11.0,tol);
         
   	ak10=500.0/tau0;

   	n10=n1+int(log(t10/taurend)/dlntau0);
         
   	if (zri!=0.0)
   	{
	   	nriend=j2ri1+nri+nri0;
      	t10=max(t10,atau0[nriend]);
      	n10=nriend+int(log(t10/atau0[nriend])/dlntau0);
      	}
      }	
	
	else
   {
	   t10=tau0;
      ak10=akmax;
      n10=nstep+1;
      }
   
   nk=0;																				// It is set to zero because it will be required later if only transfer functions are requested only.

//	FILE *fdkk;
//	fdkk = fopen("dk1.txt","w");

	// (P10) : If Cl is requested
	  
  	if (ict!=1)																		// Calculation of the CMB sources.
	{
		/*************************************************************************************************************************/
		/* set k values for which the sources for the anisotropy and polarization will be calculated. 									 */	
		/* For low values of k we use a logarithmic spacing.																							 */
		/*************************************************************************************************************************/
		
		if (zri!=0.0)
			dlnk0=2.0*dlnk;
      else
      	dlnk0=5.0*dlnk;														// if no reionization then take larger steps in the calculation			
		
		dkn1=0.8/taurmax;
      dkn2=1.5/taurmax;
			
      nk1=int(log(dkn1/(akmin*dlnk0))/dlnk0)+1;
		
		akchange=5.0*3.14159265/shor;
		
      if (akmax>akchange)
      {
      	nk2=int((akchange-akmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1;
         nk=int((akmax-akchange)/dkn2)+nk2+1;
         }
      else
      {
	      nk=int((akmax-akmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1;
         nk2=nk+1;
         }
		
      if (nk>nk0)
      {
	      printf("Sorry, the arrays were dimensioned for a maximum of");
         printf("%d, k modes. The model you requested needs, %d",nk0,nk);
         printf("Please make the arrays bigger by making ");
         printf("nk0 bigger where it appears");
         exit(1);
         }
         
		// (P11) : store k values in ak0		
      for(int ik=1;ik<=nk;ik++)
      {
	      if (ik<=nk1)
         	ak0[ik]=akmin*exp((ik-1)*dlnk0);
         else
         {
	         if(ik>nk2)
            	ak0[ik]=ak0[nk2]+(ik-nk2)*dkn2;
            else
            	ak0[ik]=ak0[nk1]+(ik-nk1)*dkn1;
            }
 			}
		
      kkmin=1;
      kkmax=nk;
		
		/************************************************************************************************************************/
		/*    K split																																				*/
		/************************************************************************************************************************/
      if((kcutflag)==1)
      {
	      for(int k=1;k<=nk;k++)
         	if (ak0[k]<(1.3*aksplit))
            	kkmax=k;
         kkmax=minint(nk,kkmax+4);
         }
         
      if ((kcutflag)==-1)
      {
	      for(int k=1;k<=nk;k++)
         {
	         if (ak0[k]<(0.3*aksplit))
            	kkmin=k;
            }   
         kkmin=maxint(1,kkmin-4);
         }
		
		/*********************************************************************************************************************************/
		/*  Main calculation. Loop over wavenumbers.																													*/
		/*  For calculation of Cl we have to calculate the source terms in the different time and different wavenumber							*/
		/*********************************************************************************************************************************/
L002:
		for(int iii=1;iii<=nvar;iii++)
			yprime[iii]=0.0;
			
		for(int ik=kkmin;ik<=kkmax;ik++)																		// Minimum and maximum value of the webnumber
		{
			ak=ak0[ik];																								// Discretized values of k's.  
																														// Begin when wave is far outside horizon.
			taustart=0.001/ak;																					// Conformal time (in Mpc) in the radiation era.
																														// Start calculations early in the radiation era. 
			taustart=min(taustart,0.1);																		// We are not including neutrinos as sources of the tensor modes.
            
   		taustart0=taustart;																					// Starting time should be saved as required for further calculations
   		
			if(amnu!=0.0)																							// Start when massive neutrinos are strongly relativistic.
   		{
	   		arel=1-3/amnu;
	   		taustart=min(taustart0,arel/adotrad);														// In radiation era a ~ t	
	   		}
		   ak2=ak*ak;																								// k^2 = k*k;
		   if(itflag!=2)
		   {																											// itflag = 2 -> Tensor perturbations are not requested  	
   			finitial(y,taustart);																			// As Cl is required so initilize the terms for calculating the source        
				
  				if(ndyn==1 || ndyn==2)																			// ndyn =1,2 -> Dark energy is perturbated
		   	{
				   y[nvar-1] = y[nvar0-1];
		     	 	y[nvar]   = y[nvar0];
		      	}
		   	tau=taustart;																						// d[] -> the sources for the anisotropy, 
																														// dp[] -> sources for the polarization,   t stands for tensor.
		   	if(ict!=0) itf=1;																					//	Transfer function is required  [transfer function loov var  1 <= itf <= nft ]
		   	if(lensflag!=0) itau=1;																			// Lensing is required	

		   	tol1=min(tol,(1.0e-17)*(ak2/1.0e-10));														// Tollerance value. The instability can ocure at low k. So increase the tolerance at low k 									
				
		   	if (ict==0)																							// No transfer functions are requested 
   				taulasttf=0.0;																						
  				else
		   	   taulasttf=tautf[ntf];																		//  It stores the time time where the last transfer function is requested	
				
				/******************************************************************************************************************************/
				/*  Integrator Communication Variables																														*/
				/******************************************************************************************************************************/
				
		   	ind=2;																								// This is specified for the "dverk" intigrator
				for(int j=1;j<=24;j++)
					c[j]=0.0;																						// This is specified for the "dverk" intigrator
	   		
			   c[3]=1.0e-8;																						// This is specified for the "dverk" intigrator
				
				nw = nvar0;
				/******************************************************************************************************************************/
				/*  Loop over timesteps.																																		*/
				/*	 This part only calculate the Scaler Source Terms and the Transfer Function																*/
				/******************************************************************************************************************************/
				for(int j=2;j<=nstep;j++)
				{
					tauend=atau0[j];	
		      	if (ak>ak10 && tauend>t10 && (ict==0 || tau>tautf[ntf]))								// ak10 & t10 srtores maximum values of k and t to be intigrated
		      	{																										// So if k or t crosses this limit then just store 0 in the source terms 		
			   	   d[ik][j]=0.0;
		      	   dp[ik][j]=0.0;
		      	   }
				    		  
    		   	else
      	 		{
	   		   	
	   		   	{
							tau = dverk(j,nvar,fderivs,tau,y,tauend,tol1,ind,c,nw,w);									// Intigrator
	   		   		ind =3;
	   		   		}
	      	      
	      	      fderivs(nvar,tau,y,yprime);
         			foutput(nvar,y,yprime,j,tau0,tau,zsanpass);
  		   	         			
         	      d[ik][j]=zsanpass[0];
         	      dp[ik][j]=zsanpass[1];
         	      dkk[ik][j]=zsanpass[2];
         	      phi=zsanpass[3];
				
						/************************************************************************************************************************/
						/*  Calculate the transfer function.																												*/
						/************************************************************************************************************************/
   		  		   do //for(;(j<nstep && itf<=ntf && atau0[j+1]>tautf[itf]);)
						{
							if (ict!=0 && itf<=ntf)						// IF NUMBER J							// itc =0 -> Only Power spectrum, ntf -> number of transfer function required
     				    	{
							
	     				    	if (j<nstep && tauend<tautf[itf] && atau0[j+1]>tautf[itf])
			   		   	{
				   		  
									tau = dverk(j,nvar,fderivs,tau,y,tautf[itf],tol1,ind,c,nw,w);		// Intigrator
	   				   		ind =3;
	   		   				}
																															// output transfer functions for this k-value.
            				if (fabs(tau-tautf[itf])<1.0e-5)		// IF NUMBER K
            				{
			   	         	if(indtau[itf]!=0)
	   	   		      		outtransf(nvar,y,0.0,indtau[itf]);										// It will calculate and store the transfer function
									else
            	   			{
            	   				if(lensflag!=0)
            	   				{
	         	      				fderivs(nvar,tau,y,yprime);
            	 	      			foutput(nvar,y,yprime,j,tau0,tau,zsanpass);						// For lensing we only need the lensing parameter phi. so store the sources in temporary variable
            	         			phi=zsanpass[3];
            	         			phik[ik][itau]=phi;
											itau=itau+1;
            	         			}						// IF END
            	      			}							// ELSE END
            	      			
	        						itf=itf+1;
	         	      		}								// IF NUMBER K END
   	      	 			}									// IF NUMBER J END	     
							}while((j<nstep) && (itf<=ntf) && (atau0[j+1]>tautf[itf]) && (ict!=0));	// FOR ENDS HERE	  
						}	                              // else end here
//					fprintf(fdkk,"%e\t",dkk[ik][j]);						
					}												// FOR OVER J WILL END HERE	
//				fprintf(fdkk,"\n");
				d[ik][nstep]=0.0;
	   		dp[ik][nstep]=0.0;
	   		}													// IF OVER ITFLAG ENDS HERE
				   																
   		/******************************************************************************************************************************/
			/*  Loop over timesteps.																																		*/
			/*	 Tensors will only be calculated if k*tau< stsp.																									*/
			/******************************************************************************************************************************/ 

   		stpt=50.0;

			if(itflag!=0)
			{
				finitialt(taustart0,yt);
   			tau=taustart0;
				
				/******************************************************************************************************************************/
				/*  Integrator Communication Variables																														*/
				/******************************************************************************************************************************/

			   ind=1;																									// This is specified for the "dverk" intigrator
				nwt=nvart+5;
			  	/***************************************************************************************************************************/
				/*  Loop over timesteps.																																	*/
				/*	 dt temp, dte electric field, dtb magnetic field		 																						*/
				/***************************************************************************************************************************/ 
   			fderivst(nvart,tau,yt,ytprime);
   			for(int j=2;j<=nstep;j++)
   			{
	   			tauend=atau0[j];
      			if((ak*tauend)>stpt)
      			{
	     				dt[ik][j]=0.0;
        				dte[ik][j]=0.0;
         			dtb[ik][j]=0.0;
         			}
     	 			else
      			{
	   		   	tol1=tol;
	   		   	{
							tau = dverk(j,nvart,fderivst,tau,yt,tauend,tol1,ind,c,nwt,wt);									// Intigrator
	   		   		ind =3;
	   		   		}
	   		   	fderivst(nvart,tau,yt,ytprime);														// Tensor derivatives
	      			foutputt(nvart,yt,ytprime,j,tau0,tau,ysanpass);									// Tensor source perturbation
	     			 	dt[ik][j] = ysanpass[0];
	     			 	dte[ik][j] = ysanpass[1];
	     			 	dtb[ik][j] = ysanpass[2];	
	     			 	}
	   			}
				dt[ik][nstep]=0.0;
		   	dte[ik][nstep]=0.0;
		   	dtb[ik][nstep]=0.0;
				}
			}						// FOR LOOP OVWE WAVE NUMBER ENDS HERE


		/*****************************************************************************************************************************/
		//	This part is added for the two field inflation models 																						  //
		/*****************************************************************************************************************************/
		
		if(twofield == 1)
		{	
			for(int storei=0;storei<nk0+1;storei++)
				for(int storej=0;storej<nstep0+1;storej++)
				{
					if(initfl == 1)
					{
						dadb[storei][storej] = d[storei][storej];
						dpadb[storei][storej] = dp[storei][storej];
						dkkadb[storei][storej] = dkk[storei][storej];
						}
					else if(initfl == 2)	
					{
						diso[storei][storej] = d[storei][storej];
						dpiso[storei][storej] = dp[storei][storej];	
						dkkiso[storei][storej] = dkk[storei][storej];
						}
					else
					{
						printf("\nERROR: In two field inflation only CDM isocurvature models are expected.");
						exit(1);
						}
					}
			if(initfl == 1)
			{
				initfl = 2;
				goto L002; 																			//	Redo the entire thing for CDM Isocurvature modes	
				}
			}
		
		/*****************************************************************************************************************************/
		}							// IF ict !=1 ENDS HERE
	


	




	printf("Integration over the perturbation variables for calculating Cl are finished");
	/******************************************************************************************************************************/
	/*  If the transfer functions are requested then ict!=0. So the transfer functions are calculated here	  							*/
	/******************************************************************************************************************************/
	int j;

	if(ict!=0)																						// If transferfunctions are requested 	 
	{
		if (ict==2) 
		{
			if (akmaxt>akmax) 
			{
				nkt=nk+int((log(akmaxt)-log(akmax))*nlnkt)+1;
     		   akdone=ak0[nk];
   	      }
   	   else
     			nkt=nk;
   		}
      	
   	else
  		{
			nkt=int((log(akmaxt)-log(akmin))*nlnkt)+1;
     		akdone=akmin;
			}
		
		akgf=1.0;
		
   	if(akdone > akgf) 
   		akgf=akdone*exp(1.5/nlnkt);								
		
		/*****************************************************************************************************************************/
		/*	 Loop over wavenumbers.																																	  */
		/*****************************************************************************************************************************/
	
		double exptemp;
   	for(int ik=nk+1;ik<=nkt;ik++)
   	{
	   	exptemp = (double)(ik-nk)/nlnkt;
	   	ak=akdone*exp(exptemp);																	// Calculate k, steps are exponentially seperated.   
  			//printf("\n exp = %e",exp(exptemp));
      		 																										// |...................|...........|
      																												//				      done(nk) 		(nkt)	
     		if(lensflag!=0)
      		itau=1;

			/******************************************************************************************************************************/
			/*  Integrator Communication Variables																														*/
			/******************************************************************************************************************************/
			
			ind =2;																									
   
         for(int jt=1;j<=24;j++)
				c[jt]=0.0;

         c[3]=1.0e-10;
				
			/******************************************************************************************************************************/
			/*	 Time loop																																						*/
			/*  Begin when wave is far outside horizon.																												*/
			/*	 Conformal time (in Mpc) in the radiation era, for photons plus 3 species of relativistic neutrinos.								*/
			/******************************************************************************************************************************/ 


      	taustart=0.001/ak;																				
      
     		taustart=min(taustart,0.1);																	// Start early in the radiation era.
      
      	if (amnu!=0.0)
      	{																										// Start when massive neutrinos are strongly relativistic.
				arel=1.0-3/amnu;										
         	taustart=min(taustart,arel/adotrad);
      		}
      		
      	ak2=ak*ak;																							// ak2 = k^2
      	finitial(y,taustart);																			// initiallize source terms in y at the time taustart
      	
      	tau=taustart;
      	
      	j = nk;

     		
 			if(ak<akgf)
			{
      		for(int itf=1;itf<=ntf;itf++)
      		{
         		if (fabs(tautf[itf]-tau)>1.0e-5)
	   		  	{
						tau = dverk(j,nvar,fderivs,tau,y,tautf[itf],tol,ind,c,nw,w);									// Intigrator
	   		   	ind =3;
	   		   	
	   		   	}
         	   else
           			tau=tautf[itf];
           				
					if(indtau[itf]!=0)																		// If transferfunction requied at this point then calculate that.
         	   	outtransf(nvar,y,0.0,indtau[itf]);
         	   else
           		{
        		    	if(lensflag!=0)																		// If lensing required then store phi	
           			{
              			if(itf==1)
              			{ 
                 			fderivs(nvar,tau,y,yprime);
                 		   foutput(nvar,y,yprime,j,tau0,tau,zsanpass);
                 		   phi1=zsanpass[3];
                    		}
							                  
 		               fderivs(nvar,tau,y,yprime);
      	            foutput(nvar,y,yprime,j,tau0,tau,zsanpass);
      	            phi=zsanpass[3];
                 		
                		
		   				phik[ik][itau]=phi;
                 		itau=itau+1;
              			}
           			}
					}
			
				closestknum=istore[1]-1;
				}
 		
 			else																									//  If k > final k specified then just store the powers and for Cl and TF
			{

	   	  	{
					tau = dverk(j,nvar,fderivs,tau,y,tauend,tol,ind,c,nw,w);									// Intigrator
	   	   	ind =3;
	   	   	}

				if(indtau[itf]!=0)																			//  If transfer function is requested at this point then call outtrans()
        			outtransf(nvar,y,0.0,indtau[itf]);

        		for(int itf=2;itf<=ntf;itf++)
        		{
	        		if(indtau[itf]!=0)
	        		{
	      		   gftemp1 = omegac*power[6*closestknum+2][indtau[itf]];						//
	      		   gftemp2 = omegab*power[6*closestknum+3][indtau[itf]];						//
	        			gftemp3 = omegan*power[6*closestknum+6][indtau[itf]];						//	
	         		gftemp4 = omegac*power[6*closestknum+2][1];									//		
	         		gftemp5 = omegab*power[6*closestknum+3][1];									//	  Growth function calculation.				
	         		gftemp6 = omegan*power[6*closestknum+6][1];									//	  Dodelson P.183  
						  																							//
	         		growthfactor=(gftemp1+gftemp2+gftemp3)/(gftemp4+gftemp5+gftemp6);		//
						
						
        		      power[6*istore[indtau[itf]+1]][itf]=power[6*istore[indtau[itf]]+1][1];								//  Power				
																																					//  All followings are transfer functions
          		   power[6*istore[itf]+2][indtau[itf]]=power[6*istore[indtau[itf]]+2][1]*growthfactor;				//
           		   power[6*istore[indtau[itf]]+3][indtau[itf]]=power[6*istore[indtau[itf]]+3][1]*growthfactor;	//
        			   power[6*istore[indtau[itf]]+4][indtau[itf]]=power[6*istore[indtau[itf]]+4][1]*growthfactor;	//
   	        		power[6*istore[indtau[itf]]+5][indtau[itf]]=power[6*istore[indtau[itf]]+5][1]*growthfactor;	//
              		power[6*istore[indtau[itf]]+6][indtau[itf]]=power[6*istore[indtau[itf]]+6][1]*growthfactor;	//
              		istore[indtau[itf]] = istore[indtau[itf]] +1;
              		}    
      	      else
           		{
	           		if (lensflag!=0)
	           		{
		           		fderivs(nvar,tau,y,yprime);
                 		foutput(nvar,y,yprime,j,tau0,tau,zsanpass);
							phihk = zsanpass[3];
                 		phik[ik][itau]=phihk*phi/phi1;
                 		itau=itau+1;
                 		}				// END IF
              		}					// END ELSE	
					}						// END FOR	
				}               		// END ELSE
			}								// END FOR WAVE NUMBER LOOP
		} 									// IF OF ICT!=0 ENDS HERE      
	

	if (lensflag!=0) 
		ntf=ntf-ntau;

	if(ict!=0)
		output_power(ntf,amnu);
	

//fclose(fdkk);

	// (P13) : Calculating the Cl
	
	/*********************************************************************************************************************************/
	/*  If CMB calculations are requested, calculate the Cl by integrating the sources over time and over k. 								*/
	/*********************************************************************************************************************************/
			
	if (ict != 1)
	{
		/******************************************************************************************************************************/
		/*  If Scaler Cl's are requested then initiallize the Cls and calculate the derivative of the source terms.							*/
		/******************************************************************************************************************************/
		printf("l0 = %d",l0);
			
		if (itflag!=2)
		{
			if(twofield == 0)
			{
				for(int i=1;i<=nstep;i++)
			   {
				   for(int ii=0;ii < nk0+1 ;ii++)
				   {
					   sand[ii]=d[ii][i];
					   sandp[ii]=dp[ii][i];
					   sandkk[ii]=dkk[ii][i];
					   }
			   	spline(ak0,sand,nk,d0lo,d0hi,sandpr);											//									
  		      	spline(ak0,sandp,nk,d0lo,d0hi,sandppr);										// Calculate the space derivative of the 
     		   	spline(ak0,sandkk,nk,d0lo,d0hi,sandkpr);										// scaler source terms
     		   	for(int ii=0;ii < nk0+1 ;ii++)
     		   	{
	     			   dpr[ii][i]=sandpr[ii];
	     			   dppr[ii][i]=sandppr[ii];
	     			   dkpr[ii][i]=sandkpr[ii];
	     			   }
     		   	}
				}
			
			/***************************************************************************************************************************************/
			// This part is added for two field inflation 																														//
			/***************************************************************************************************************************************/

			else if(twofield == 1)																		// For two field model
			{
				for(int i=1;i<=nstep;i++)
			   {
				   for(int ii=0;ii < nk0+1 ;ii++)
				   {
					   sand[ii]=dadb[ii][i];
					   sandp[ii]=dpadb[ii][i];
					   sandkk[ii]=dkkadb[ii][i];
					   }
			   	spline(ak0,sand,nk,d0lo,d0hi,sandpr);											//									
  		      	spline(ak0,sandp,nk,d0lo,d0hi,sandppr);										// Calculate the space derivative of the 
     		   	spline(ak0,sandkk,nk,d0lo,d0hi,sandkpr);										// scaler source terms
     		   	for(int ii=0;ii < nk0+1 ;ii++)
     		   	{
	     			   dpradb[ii][i]=sandpr[ii];
	     			   dppradb[ii][i]=sandppr[ii];
	     			   dkpradb[ii][i]=sandkpr[ii];
	     			   }
     		   	}
				for(int i=1;i<=nstep;i++)
			   {
				   for(int ii=0;ii < nk0+1 ;ii++)
				   {
					   sand[ii]=diso[ii][i];
					   sandp[ii]=dpiso[ii][i];
					   sandkk[ii]=dkkiso[ii][i];
					   }
			   	spline(ak0,sand,nk,d0lo,d0hi,sandpr);											//									
  		      	spline(ak0,sandp,nk,d0lo,d0hi,sandppr);										// Calculate the space derivative of the 
     		   	spline(ak0,sandkk,nk,d0lo,d0hi,sandkpr);										// scaler source terms
     		   	for(int ii=0;ii < nk0+1 ;ii++)
     		   	{
	     			   dpriso[ii][i]=sandpr[ii];
	     			   dppriso[ii][i]=sandppr[ii];
	     			   dkpriso[ii][i]=sandkpr[ii];
	     			   }
     		   	}
				}
			/***************************************************************************************************************************************/



//			for(int in=1;in<=nn;in++)
//			{
//				for(int j=1;j<=l0;j++)
//				{
//											
//					cl[j][in]=0.0;																				//	
//      		   cpl[j][in]=0.0;																			// Initiallize the power spectrum terms
//   			   ccl[j][in]=0.0;																			// 
//					ckkl[j][in]=0.0;																			//
//         		ctkl[j][in]=0.0;																			//		
//           		}
//        		}
     		}
		/******************************************************************************************************************************/
		/*  If Tensor Cl's are requested then initiallize the Cls and calculate the derivative of the source terms.							*/
		/******************************************************************************************************************************/
		
	   if (itflag!=0)
  		{
   		for(int i=1;i<=nstep;i++)
  			{
			   for(int ii=0;ii < nk0+1 ;ii++)
			   {
				   sandt[ii]=dt[ii][i];
				   sandte[ii]=dte[ii][i];
				   sandtb[ii]=dtb[ii][i];
//					printf("%e\t%e\t%e \n",sandt[ii],sandte[ii],sandtb[ii]);
				   }
//				exit(1);   
				spline(ak0,sandt,nk,d0lo,d0hi,sandtpr);										//
        		spline(ak0,sandte,nk,d0lo,d0hi,sandtepr);										// Calculate the time derivative of the  
        		spline(ak0,sandtb,nk,d0lo,d0hi,sandtbpr);										// tensor source terms
     		   for(int ii=0;ii < nk0+1 ;ii++)
     		   {
	     		   dtpr[ii][i]=sandtpr[ii];
	     		   dtepr[ii][i]=sandtepr[ii];
	     		   dtbpr[ii][i]=sandtbpr[ii];
	     		   }
        		}
        	
//   		for(int in=1;in<=nn;in++)
//   		{
//   		   for(int j=1;j<=l0;j++)
//      		{
//	      		ctl[j][in]=0.0;																			//
//         		ctel[j][in]=0.0;																			// Initiallize the tensor power spectrum
//         		ctbl[j][in]=0.0;																			//
//         		ctcl[j][in]=0.0;																			//	
//         		}
//      		}
     		}   

		// (P14) : Define the wavenumber and the difference drid for integration      
		/******************************************************************************************************************************/
		/*  Fixing the wavenumbers for intigration.																												*/
		/******************************************************************************************************************************/         
		
		dk=2.5/tau0;																							
  		dk0=1.5/tau0;
  		no=700;																									// The dimension of K arrays rellated to no =700		
  		dlnk1=0.07;
		
  		no1=(int)(log(10.0*dk0/akmin)/dlnk1)+1;
 		
  		if (akmax > (no*dk0))																				
  			nko=(int)((akmax-no*dk0)/dk)+no;
  		else
  		{
			no=(int)((akmax-10.0*dk0)/dk0)+no1;
     		nko=no;
     		}
         
		if(nko > nkmax)																						// Dimension of K array is nkmax
  		{
			printf("\nSorry, the arrays were dimensioned for a maximum of %d k modes", nkmax); 
     		printf("\nThe model you requested needs %d",nko);
     		printf("\nPlease make the arrays bigger by making nkmax bigger where it appears");
     		exit(1);
     		}
			
		
		for(int k=1;k<=nko;k++)																				//   |...................|...............|..........................|
  		{																											//		 <for short wev> (no1) <mid wave> (no)  <long wave length>
		   if (k<=no)																							//        log space            linear         linear large space 
   		{																										//	
	   		if (k<=no1)																						//    Wave number grid points are stored in ak[]	
	  		 	{																									//    Wave number spacing are stoder in dak1[]	
		   		ak1[k]=10.0*dk0*exp(-(no1-k)*dlnk1);												//
           		dak1[k]=ak1[k]*dlnk1;																	//
           		}
               
        		else
        		{
         		ak1[k]=ak1[no1]+(k-no1)*dk0;
           		dak1[k]=dk0;
           		}
        		}
                       
     		else
     		{
      		ak1[k]=ak1[no]+(k-no)*dk;  	
        		dak1[k]=dk;
        		}
     		}   
					
		dak1[1]=0.5*dak1[1];
  		dak1[no1]=0.5*(dak1[no1]+dk0);
  		dak1[nko]=0.5*dak1[nko];

		/*****************************************************************************************************************************/
		/*   K split																																					  */
		/*****************************************************************************************************************************/
			
																														// |Kmin..............................|Kmax
																															
		kkmin=1;																										// Minimum Value of K (Wavenumber)			
		kkmax=nko;																									// Maximum value of K (Wavenumber)	
  		 
  	  		    
		if (kcutflag==1)																							// |Kmin............|Kmax.................|
		{																												// |1...............|1.3..................|							
			for(int k=1;k<=nko;k++)																				// K < K*
				if(ak1[k]<1.3*aksplit)kkmax=k+1;															
   		}
      
		if (kcutflag==-1)																							// |............|Kmin.................|Kmax
		{																												// |............|0.3..................|nko	
			for(int k=1;k<=nko;k++)																				// K > K*
				if(ak1[k]<=0.3*aksplit)kkmin=k;
			}

		/*****************************************************************************************************************************/
		/*   Start the loops on wavenumber. This loop will interpolate source terms for scaler Cl	at the discritized nodes			  */
		/*****************************************************************************************************************************/

		klo=1;
		khi=2;
		
						
		for(int k=kkmin;k<=kkmax;k++)
		{
			akt=ak1[k];																								//	Finding position of k in table ak0 to do the interpolation.
  			
			for(klo=1,khi=2; ((akt > ak0[klo+1]) && (klo < nk-1)); klo++);							// Check.. There may be some problem.
  			khi=klo+1;
			
			
			ho=ak0[khi]-ak0[klo];																				
		 	a0=(ak0[khi]-akt)/ho;
			b0=(akt-ak0[klo])/ho;
			
	  		nstps=0;
  			nstpt=0;
			
			/***********************************************************************************************************************/
			/*     If scalers wanted then this loop will interpolate the source terms at the nodes for tensor Cl						  */	
			/***********************************************************************************************************************/
			
  			if (itflag!=2)																							// itflag=2 -> Only Tensor part
  			{   
   			for(int j=1;j<=l0;j++)
   			{
					if(twofield == 0)
					{
		   			dl[j]=0.0;
   	     			dl2[j]=0.0;
   	     			dl3[j]=0.0;
  				    	dpl[j]=0.0;
   	  			 	dpl2[j]=0.0;
   	  			 	dpl3[j]=0.0;
   	    			dkl[j]=0.0;
   	    			dkl2[j]=0.0;
   	     			dkl3[j]=0.0;
						}

					else if(twofield == 1)
					{
		   			dladb[j]=0.0;
   	     			dl2adb[j]=0.0;
   	     			dl3adb[j]=0.0;
  				    	dpladb[j]=0.0;
   	  			 	dpl2adb[j]=0.0;
   	  			 	dpl3adb[j]=0.0;
   	    			dkladb[j]=0.0;
   	    			dkl2adb[j]=0.0;
   	     			dkl3adb[j]=0.0;

		   			dliso[j]=0.0;
   	     			dl2iso[j]=0.0;
   	     			dl3iso[j]=0.0;
  				    	dpliso[j]=0.0;
   	  			 	dpl2iso[j]=0.0;
   	  			 	dpl3iso[j]=0.0;
   	    			dkliso[j]=0.0;
   	    			dkl2iso[j]=0.0;
   	     			dkl3iso[j]=0.0;
						}
					}
	
																														// Interpolating the source as a function of time for the present wavelength.
     			if(akt<ak10)
     				nstps=nstep-1;
     			else
        			nstps=n10;
       		       
				if(twofield == 0)
				{
					s2[1]=0.0;																							// Interpolating the source terms	
   	  			sp2[1]=0.0;																							// The interpolated source terms will be temporarily stored in these variables		
   	  			sk2[1]=0.0;																							//
     				}
		
				else if(twofield == 1)
				{
					s2adb[1]=0.0;																							// Interpolating the source terms	
   	  			sp2adb[1]=0.0;																							// The interpolated source terms will be temporarily stored in these variables		
   	  			sk2adb[1]=0.0;																							//

					s2iso[1]=0.0;																							// Interpolating the source terms	
   	  			sp2iso[1]=0.0;																							// The interpolated source terms will be temporarily stored in these variables		
   	  			sk2iso[1]=0.0;																							//
					}
     			
  			 	for(int i=2;i<=nstps;i++)																				// Spline interpolation of the source terms 
     			{
					if(twofield == 0)
					{
	      			s2[i]=a0*d[klo][i]+b0*d[khi][i]+((a0*a0*a0-a0)*dpr[klo][i]+(b0*b0*b0-b0)*dpr[khi][i])*ho*ho/6.0;
	        			sp2[i]=a0*dp[klo][i]+b0*dp[khi][i]+((a0*a0*a0-a0)*dppr[klo][i]+(b0*b0*b0-b0)*dppr[khi][i])*ho*ho/6.0;
	        			sk2[i]=a0*dkk[klo][i]+b0*dkk[khi][i]+((a0*a0*a0-a0)*dkpr[klo][i]+(b0*b0*b0-b0)*dkpr[khi][i])*ho*ho/6.0;
						}

					else if(twofield == 1)
					{
	      			s2adb[i]=a0*dadb[klo][i]+b0*dadb[khi][i]+((a0*a0*a0-a0)*dpradb[klo][i]+(b0*b0*b0-b0)*dpradb[khi][i])*ho*ho/6.0;
	        			sp2adb[i]=a0*dpadb[klo][i]+b0*dpadb[khi][i]+((a0*a0*a0-a0)*dppradb[klo][i]+(b0*b0*b0-b0)*dppradb[khi][i])*ho*ho/6.0;
	        			sk2adb[i]=a0*dkkadb[klo][i]+b0*dkkadb[khi][i]+((a0*a0*a0-a0)*dkpradb[klo][i]+(b0*b0*b0-b0)*dkpradb[khi][i])*ho*ho/6.0;

	      			s2iso[i]=a0*diso[klo][i]+b0*diso[khi][i]+((a0*a0*a0-a0)*dpriso[klo][i]+(b0*b0*b0-b0)*dpriso[khi][i])*ho*ho/6.0;
	        			sp2iso[i]=a0*dpiso[klo][i]+b0*dpiso[khi][i]+((a0*a0*a0-a0)*dppriso[klo][i]+(b0*b0*b0-b0)*dppriso[khi][i])*ho*ho/6.0;
	        			sk2iso[i]=a0*dkkiso[klo][i]+b0*dkkiso[khi][i]+((a0*a0*a0-a0)*dkpriso[klo][i]+(b0*b0*b0-b0)*dkpriso[khi][i])*ho*ho/6.0;
						}
					}

				if(twofield == 0)
				{
	      		s2[nstps+1]=0.0;
  		   		sp2[nstps+1]=0.0;
   	  			sk2[nstps+1]=0.0;
					}

				else if(twofield == 1)
				{
					s2adb[nstps+1]=0.0;																	
   	  			sp2adb[nstps+1]=0.0;																		
   	  			sk2adb[nstps+1]=0.0;																

					s2iso[nstps+1]=0.0;																	
   	  			sp2iso[nstps+1]=0.0;																		
   	  			sk2iso[nstps+1]=0.0;																
					}
				}
                                                                                          			
			/***********************************************************************************************************************/
			/*     If tensors wanted then this loop will interpolate the source terms at the nodes for tensor Cl						  */	
			/***********************************************************************************************************************/
			
			if (itflag!=0)
 			{
   			for(int j=1;j<=l0;j++)
   			{
					dtl[j]=0.0;
     				dtl2[j]=0.0;
     				dtl3[j]=0.0;
     				dtel2[j]=0.0;
     				dtbl2[j]=0.0;
	  				}
	  			
																													// Interpolating the tensor source as a function of time for the present wavelength.
     			st2[1]=0.0;
     			ste2[1]=0.0;
     			stb2[1]=0.0;																					// Initiallize the source terms
     			nstpt=2;
      			     
     			for(int i=2;i<=nstep;i++)
     			{
   		  		xf=akt*(tau0-atau0[i]);
     				if(((akt*atau0[i])<stpt) && (xf>1.0e-8))											// Spline interpolation of the source terms 
     				{
        				nstpt=i;
        				st2[i]=a0*dt[klo][i]+b0*dt[khi][i]+((a0*a0*a0-a0)*dtpr[klo][i]+(b0*b0*b0-b0)*dtpr[khi][i])*ho*ho/6.0;
        				ste2[i]=a0*dte[klo][i]+b0*dte[khi][i]+((a0*a0*a0-a0)*dtepr[klo][i]+(b0*b0*b0-b0)*dtepr[khi][i])*ho*ho/6.0;
        				stb2[i]=a0*dtb[klo][i]+b0*dtb[khi][i]+((a0*a0*a0-a0)*dtbpr[klo][i]+(b0*b0*b0-b0)*dtbpr[khi][i])*ho*ho/6.0;
        				}
         		        
     	 			else
     				{
        				st2[i]=0.0;
        				ste2[i]=0.0;
        				stb2[i]=0.0;
        				}	
					}
				
	      	nstpt=(nstpt>n1)?nstpt:n1;
	      	st2[nstpt]=0.0;
		   	ste2[nstpt]=0.0;
		   	stb2[nstpt]=0.0;
     			}	
			
			// (P16) : store k*tau in mxx() for calculating the bessel functions
			
     		/***********************************************************************************************************************/
			/*    Findind the position in the xx table for the x correponding to each timestep												  */	
			/***********************************************************************************************************************/
 		 	
			for(int i=1;i<=maxint(nstps,nstpt)+2;i++)
			{
				xf=fabs(akt*(tau0-atau0[i]));																// Store k*t values at the grid points	
				
     			if(xf<=5.0)
     				de2=10.0*xf+1.0;               
      		
     			else if (xf<=25.0 && xf>5.0)
     				de2=(xf-5.0)*5.0+51.0;
				
     			else
     				de2=(xf-25.0)+151.0;
				
  		   	mxx[i]=(int)(de2);
     			mxx[i]=maxint(mxx[i],1);																	// Tempurarily store k*t values at the grid points
				}
				
			/*********************************************************************************************************/				
			/*     Begin l and  time-loop to integrate scalar perturbations.														*/	
			/*     Determining ranges of integration																						*/
			/*********************************************************************************************************/						
			
     		// (P17) : Calculating the Cl
     		
     		if (itflag!=2)
     		{
   		   for(int j=1;j<=l0;j++)
        		{
         		xlim=0.05*l[j];
           		xlim=max(xlim,xlimmin);
           		xlim=l[j]-xlim;
           		xlmax1=80.0*l[j];
           		xlmax2=min(2.0*l[j],1.0*kmaxjl);
					
  		         tmin=tau0-xlmax1/akt;                      
     		      tmax=tau0-xlim/akt;
        		   tmax=min(tau0,tmax);
        		   tmin=max(atau0[2],tmin);
					
  		         if(tmax<atau0[2]) 
  		         {
  		         	goto MOVETOA;
						}
  		         if(zri==0.0)
					{
						if(tmin<taurend)
           		      nstart1=2;
           		   else
              			nstart1=n1+(int)(log(tmin/taurend)/dlntau0);
                    
           			if (tmax<taurend)
              		{
                 		nstop1=n1;
                 		nstop1a=n1;
                 		}
               		
              		else
              		{
              		   nstop1=n1+(int)(log(tmax/taurend)/dlntau0);
              		   nstop1=minint(nstop1,nstps);
              		    
              		   if ((akt*dtau1[nstop1]) > 1.0)
              				nstop1a=maxint(n1,n1-int(log(akt*dlntau0*taurend)/dlntau0));
                 		else
	              		   nstop1a=nstop1;
              			}
						}
							
            	else
            	{
	         	   if(tmin<taurend)
            	   {
              		   nstart1=2;
              		   nstart2=j2ri1;                   
              		   }
              		else
              		{
	              		if(tmin < atau0[nriend])
                 		{
	              		   nstart1=minint(n1+(int)(log((tmin/taurend))/dlntau0),j2ri1);
	              		   nstart2=j2ri1;
                 		   }
              			else
                 		{
	                 		nstart1=nstep+1;
                    		nstart2=nriend+(int)(log(tmin/atau0[nriend]/dlntau0));
                    		}
                 		}
							
              		if(tmax < taurend)
              		{
	           		   nstop1=n1;
              		   nstop2=0;
                 		}
                	
         			else
         			{  
            			if(tmax < atau0[j2ri1])
            			{
              				nstop1=n1+(int)(log(tmax/taurend)/dlntau0);
              				nstop2=0;
            				}
            			
            			else
							{
              				nstop1=j2ri1-1;
              				nstop2=maxint(nriend,nriend+(int)(log(tmax/atau0[nriend])/dlntau0));
					
              				nstop2=minint(nstop2,nstps);
                 		         
              				if((akt*dtau1[nstop2]) > 1)
              					nstop2a=maxint(nriend,nriend-(int)(log(akt*dlntau0*atau0[nriend])/dlntau0));
              				else
                 				nstop2a=nstop2;
            				}
				
            			nstop1=minint(nstop1,j2ri1-1);
                       
							if((akt*dtau1[nstop1])>1)
            				nstop1a=maxint(n1,n1-(int)(log(akt*dlntau0*taurend)/dlntau0));
            			else
            				nstop1a=nstop1;
							}
						}
															
					/***************************************************************************************************************************/
					/*    Integration before reionization																													*/
					/*		Interpolating jls at points where the sources are recorded.																				*/
					/***************************************************************************************************************************/			
					
     				for(int i=nstart1;i<=nstop1a;i++)
     				{
	     				xf=akt*(tau0-atau0[i]);
        				m2=mxx[i];                           
        				h2=xx[m2+1]-xx[m2];
        				a2=(xx[m2+1]-xf)/h2;
        				b2=(xf-xx[m2])/h2;
        				ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
        		
						if(twofield == 0)
						{	
        					dl2[j]=dl2[j]+s2[i]*ajl0*dtau2[i];
        					dpl2[j]=dpl2[j]+sp2[i]*ajl0*dtau2[i];
        					dkl2[j]=dkl2[j]+sk2[i]*ajl0*dtau2[i];
        					}

						if(twofield == 1)
						{	
        					dl2adb[j]=dl2adb[j]+s2adb[i]*ajl0*dtau2[i];
        					dpl2adb[j]=dpl2adb[j]+sp2adb[i]*ajl0*dtau2[i];
        					dkl2adb[j]=dkl2adb[j]+sk2adb[i]*ajl0*dtau2[i];

        					dl2iso[j]=dl2iso[j]+s2iso[i]*ajl0*dtau2[i];
        					dpl2iso[j]=dpl2iso[j]+sp2iso[i]*ajl0*dtau2[i];
        					dkl2iso[j]=dkl2iso[j]+sk2iso[i]*ajl0*dtau2[i];
        					}
						}
					
					/***************************************************************************************************************************/
					/*     Scaler case.																																			*/
					/*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
					/***************************************************************************************************************************/		
//					printf("\n %d  %d  ",nstop1a,nstop1); 
					if(twofield == 0)
					{
						for(int i=nstop1a+1;i<=nstop1;i++)
						{
							xf=akt*(tau0-atau0[i]);
   	      			m2=mxx[i];
   	  	    			dtau3=dtau1[i]*fabs(s2[i]/(s2[i+1]-s2[i]+1.0e-10));
//							printf("")
//			       			printf("%e %e %e %e\n",xf,xlmax2,akt,dtau3);
		        			if ((xf < xlmax2)||((akt*dtau3) < 1.0))
  		       			{
		
  	 	    	   			xi=xf-akt*dtau1[i];
     		   	 			m1=mxx[i+1];
           					ddt1=s2[i];
           					ddt2=s2[i+1];
           					ddp1=sp2[i];
           					ddp2=sp2[i+1];
           					ddk1=sk2[i];
           					ddk2=sk2[i+1];
//          	 				printf("%e %e %e\n",ddt1,ddp1,ddk1);
        	 					for(int lx=m1+1;lx<=m2;lx++)												// convert the source terms from function of (k,t) to the function of (j,l)
           					{																					// M.Zaldarriaga's thesis Eq(3.37)
           	   				x=xx[lx];
           	   				ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2;
           	   				ddp=(ddp1-ddp2)*(x-xi)/(xf-xi)+ddp2;
           	   				ddk=(ddk1-ddk2)*(x-xi)/(xf-xi)+ddk2;
           	   				dl3[j]=dl3[j]+ajl[lx][j]*ddt*dxx[lx];
          	    				dpl3[j]=dpl3[j]+ajl[lx][j]*ddp*dxx[lx];
           	   				dkl3[j]=dkl3[j]+ajl[lx][j]*ddk*dxx[lx];
           						}
        						}
        					}
						}	

					if(twofield == 1)
					{
						for(int i=nstop1a+1;i<=nstop1;i++)
						{
							xf=akt*(tau0-atau0[i]);
   	      			m2=mxx[i];
   	  	    			dtau3=dtau1[i]*fabs(s2adb[i]/(s2adb[i+1]-s2adb[i]+1.0e-10));
//							printf("")
//			       			printf("%e %e %e %e\n",xf,xlmax2,akt,dtau3);
		        			if ((xf < xlmax2)||((akt*dtau3) < 1.0))
  		       			{
		
  	 	    	   			xi=xf-akt*dtau1[i];
     		   	 			m1=mxx[i+1];

           					ddt1adb=s2adb[i];
           					ddt2adb=s2adb[i+1];
           					ddp1adb=sp2adb[i];
           					ddp2adb=sp2adb[i+1];
           					ddk1adb=sk2adb[i];
           					ddk2adb=sk2adb[i+1];

           					ddt1iso=s2iso[i];
           					ddt2iso=s2iso[i+1];
           					ddp1iso=sp2iso[i];
           					ddp2iso=sp2iso[i+1];
           					ddk1iso=sk2iso[i];
           					ddk2iso=sk2iso[i+1];
//          	 				printf("%e %e %e\n",ddt1,ddp1,ddk1);
        	 					for(int lx=m1+1;lx<=m2;lx++)												// convert the source terms from function of (k,t) to the function of (j,l)
           					{																					// M.Zaldarriaga's thesis Eq(3.37)
           	   				x=xx[lx];

           	   				ddtadb=(ddt1adb-ddt2adb)*(x-xi)/(xf-xi)+ddt2adb;
           	   				ddpadb=(ddp1adb-ddp2adb)*(x-xi)/(xf-xi)+ddp2adb;
           	   				ddkadb=(ddk1adb-ddk2adb)*(x-xi)/(xf-xi)+ddk2adb;
           	   				dl3adb[j]=dl3adb[j]+ajl[lx][j]*ddt*dxx[lx];
          	    				dpl3adb[j]=dpl3adb[j]+ajl[lx][j]*ddpadb*dxx[lx];
           	   				dkl3adb[j]=dkl3adb[j]+ajl[lx][j]*ddkadb*dxx[lx];

           	   				ddtiso=(ddt1iso-ddt2iso)*(x-xi)/(xf-xi)+ddt2iso;
           	   				ddpiso=(ddp1iso-ddp2iso)*(x-xi)/(xf-xi)+ddp2iso;
           	   				ddkiso=(ddk1iso-ddk2iso)*(x-xi)/(xf-xi)+ddk2iso;
           	   				dl3iso[j]=dl3iso[j]+ajl[lx][j]*ddtiso*dxx[lx];
          	    				dpl3iso[j]=dpl3iso[j]+ajl[lx][j]*ddpiso*dxx[lx];
           	   				dkl3iso[j]=dkl3iso[j]+ajl[lx][j]*ddkiso*dxx[lx];
           						}
        						}
        					}
						}	


					
//					printf("%d %d %e %e\n",k,j,dkl3[j],dl3[j]);         	
					
					/***************************************************************************************************************************/
					/*    Integration after reionization																													*/
					/*		Interpolating jls at points where the sources are recorded.																				*/
					/***************************************************************************************************************************/
	
     				if (zri!=0.0)
     				{
						if(twofield == 0)
						{
	     					for(int i=nstart2;i<=nstop2a;i++)
   	     				{
		        				xf=akt*(tau0-atau0[i]);
   	        				m2=mxx[i];
   	        				h2=xx[m2+1]-xx[m2];
   	        				a2=(xx[m2+1]-xf)/h2;
   	        				b2=(xf-xx[m2])/h2;
   	        				ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
   	        				dl2[j]=dl2[j]+s2[i]*ajl0*dtau2[i];
   	        				dpl2[j]=dpl2[j]+sp2[i]*ajl0*dtau2[i];
   	        				dkl2[j]=dkl2[j]+sk2[i]*ajl0*dtau2[i];
   	        				}
          				}

						if(twofield == 1)
						{
	     					for(int i=nstart2;i<=nstop2a;i++)
   	     				{
		        				xf=akt*(tau0-atau0[i]);
   	        				m2=mxx[i];
   	        				h2=xx[m2+1]-xx[m2];
   	        				a2=(xx[m2+1]-xf)/h2;
   	        				b2=(xf-xx[m2])/h2;
   	        				ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
   	        				dl2adb[j]=dl2adb[j]+s2adb[i]*ajl0*dtau2[i];
   	        				dpl2adb[j]=dpl2adb[j]+sp2adb[i]*ajl0*dtau2[i];
   	        				dkl2adb[j]=dkl2adb[j]+sk2adb[i]*ajl0*dtau2[i];
   	        				}
	     					for(int i=nstart2;i<=nstop2a;i++)
   	     				{
		        				xf=akt*(tau0-atau0[i]);
   	        				m2=mxx[i];
   	        				h2=xx[m2+1]-xx[m2];
   	        				a2=(xx[m2+1]-xf)/h2;
   	        				b2=(xf-xx[m2])/h2;
   	        				ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
   	        				dl2iso[j]=dl2iso[j]+s2iso[i]*ajl0*dtau2[i];
   	        				dpl2iso[j]=dpl2iso[j]+sp2iso[i]*ajl0*dtau2[i];
   	        				dkl2iso[j]=dkl2iso[j]+sk2iso[i]*ajl0*dtau2[i];
   	        				}
          				}

           			/***************************************************************************************************************************/
						/*     Scaler case.																																			*/
						/*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
						/***************************************************************************************************************************/		

						for(int i=nstop2a+1;i<=nstop2;i++)
        				{
	        				xi=xf-akt*dtau1[i];
           				xf=akt*(tau0-atau0[i]);
          				m2=mxx[i];
          				
							if(twofield == 0)
 	          				dtau3=dtau1[i]*fabs(s2[i]/(s2[i+1]-s2[i]+1.0e-10));

							else if(twofield == 0)															// Check may be something wrong here..			
 	          				dtau3=dtau1[i]*fabs(s2adb[i]/(s2adb[i+1]-s2adb[i]+1.0e-10));
							
           				if ((xf < xlmax2)||((akt*dtau3)<1.0))
           				{
           					m1=mxx[i+1];

              				if(twofield == 0)
								{
									ddt1=s2[i];
	              				ddt2=s2[i+1];
	              				ddp1=sp2[i];
	              				ddp2=sp2[i+1];
	              				ddk1=sk2[i];
	              				ddk2=sk2[i+1];
									}	
              				if(twofield == 0)
								{
									ddt1adb=s2adb[i];
	              				ddt2adb=s2adb[i+1];
	              				ddp1adb=sp2adb[i];
	              				ddp2adb=sp2adb[i+1];
	              				ddk1adb=sk2adb[i];
	              				ddk2adb=sk2adb[i+1];

									ddt1iso=s2iso[i];
	              				ddt2iso=s2iso[i+1];
	              				ddp1iso=sp2iso[i];
	              				ddp2iso=sp2iso[i+1];
	              				ddk1iso=sk2iso[i];
	              				ddk2iso=sk2iso[i+1];
									}

              				for(int lx=m1+1;lx<=m2;lx++)												// convert the source terms from function of (k,t) to the function of (j,l)
              				{																					// M.Zaldarriaga's thesis Eq(3.37)
            					x=xx[lx];
									if(twofield == 0)
									{
                    				ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2;
                    				ddp=(ddp1-ddp2)*(x-xi)/(xf-xi)+ddp2;
                    				ddk=(ddk1-ddk2)*(x-xi)/(xf-xi)+ddk2;
                    				dl3[j]=dl3[j]+ajl[lx][j]*ddt*dxx[lx];
                    				dpl3[j]=dpl3[j]+ajl[lx][j]*ddp*dxx[lx];
                    				dkl3[j]=dkl3[j]+ajl[lx][j]*ddk*dxx[lx];
										}	 																						
									if(twofield == 1)
									{
                    				ddtadb=(ddt1adb-ddt2adb)*(x-xi)/(xf-xi)+ddt2adb;
                    				ddpadb=(ddp1adb-ddp2adb)*(x-xi)/(xf-xi)+ddp2adb;
                    				ddkadb=(ddk1adb-ddk2adb)*(x-xi)/(xf-xi)+ddk2adb;
                    				dl3adb[j]=dl3adb[j]+ajl[lx][j]*ddtadb*dxx[lx];
                    				dpl3adb[j]=dpl3adb[j]+ajl[lx][j]*ddpadb*dxx[lx];
                    				dkl3adb[j]=dkl3adb[j]+ajl[lx][j]*ddkadb*dxx[lx];

                    				ddtiso=(ddt1iso-ddt2iso)*(x-xi)/(xf-xi)+ddt2iso;
                    				ddpiso=(ddp1iso-ddp2iso)*(x-xi)/(xf-xi)+ddp2iso;
                    				ddkiso=(ddk1iso-ddk2iso)*(x-xi)/(xf-xi)+ddk2iso;
                    				dl3iso[j]=dl3iso[j]+ajl[lx][j]*ddtiso*dxx[lx];
                    				dpl3iso[j]=dpl3iso[j]+ajl[lx][j]*ddpiso*dxx[lx];
                    				dkl3iso[j]=dkl3iso[j]+ajl[lx][j]*ddkiso*dxx[lx];
										}
             					}						// LOOP OVER lx ENDS HERE
            				}							// IF END HERE
							}								// FOR LOOP OVER I ENDS HERE 
						}									// IF OF ZRI ENDS HERE
					

					if(twofield == 0)
					{
		     			dl[j]=dl2[j]+dl3[j]/akt;
	      			dpl[j]=dpl2[j]+dpl3[j]/akt;
						dkl[j]=dkl2[j]+dkl3[j]/akt;
				//		printf("Hi, %e\n",dl[j]);
						}
					if(twofield == 1)
					{
		     			dladb[j]=dl2adb[j]+dl3adb[j]/akt;
	      			dpladb[j]=dpl2adb[j]+dpl3adb[j]/akt;
						dkladb[j]=dkl2adb[j]+dkl3adb[j]/akt;

		     			dliso[j]=dl2iso[j]+dl3iso[j]/akt;
	      			dpliso[j]=dpl2iso[j]+dpl3iso[j]/akt;
						dkliso[j]=dkl2iso[j]+dkl3iso[j]/akt;
						//printf("%e %e\n",dladb[j],dliso[j]);
						}

					//	printf("%d\t %d\t %e\t %e\t %e\n",k,j,dkl2[j],dkl3[j],dkl[j]);
					}	
														//FOR LOOP OVER J ENDS HERE
MOVETOA:

				/***************************************************************************************************************************/
				/*     Scaler case.																																			*/
				/*		 Finally calculate Cl by intigrating over K's . ps: l*(l-1) will be done in the next step										*/
				/***************************************************************************************************************************/		
				
				//printf("nn = %d %d %d\n",nn,l0,twofield);
   			for(int in=1;in<=nn;in++)
				{
					//printf("Hi.. \n");
					if(twofield == 0)
					{
		     			for(int j=1;j<=l0;j++)
						{
						//	printf("\tHi\n");
	   	   			powersflat(akt,in,&apowers);																	// apower is the initial power spectrum
         	  			apowers=apowers/akt;
           			
           				ckj=apowers*dl[j]*dl[j]*dak1[k];
           				cl[j][in]=cl[j][in]+ckj;
           			
        					cpkj=apowers*dpl[j]*dpl[j]*dak1[k];															// Calculating Cl's 	
           				cpl[j][in]=cpl[j][in]+cpkj;																	// Zaldarriaga PhD Thesis 3.25	
            		
           				cckj=apowers*dl[j]*dpl[j]*dak1[k];
        				 	ccl[j][in]=ccl[j][in]+cckj;
            		
           				ckkkj=apowers*dkl[j]*dkl[j]*dak1[k];
           				ckkl[j][in]=ckkl[j][in]+ckkkj;
           	 		
        	 				ctkkj=apowers*dl[j]*dkl[j]*dak1[k];
           				ctkl[j][in]=ctkl[j][in]+ctkkj;
							//printf("%e\t%e\t%e\%e\n",cl[j][in],apowers,dl[j],dak1[k]);
		     				}
						}

					else if(twofield == 1)
					{
						if(powerfileflag == 1)
							interpoletepowers(akt);
						//printf("Hi %d\n",twofield);
		     			for(int j=1;j<=l0;j++)
						{
	   	   			//powersflat(akt,in,&apowers);
							if(powerfileflag == 0)
								apowers = poweradiabatic(akt);																	// apower is the initial power spectrum
         	  			else if(powerfileflag == 1)
								apowers = powerA;
		
							apowers=apowers/akt;
           			
           				ckj=apowers*dladb[j]*dladb[j]*dak1[k];
           				cladb[j][in]=cladb[j][in]+ckj;
           			
        					cpkj=apowers*dpladb[j]*dpladb[j]*dak1[k];															// Calculating Cl's 	
           				cpladb[j][in]=cpladb[j][in]+cpkj;																	// Zaldarriaga PhD Thesis 3.25	
            		
           				cckj=apowers*dladb[j]*dpladb[j]*dak1[k];
        				 	ccladb[j][in]=ccladb[j][in]+cckj;
		     				}

		     			for(int j=1;j<=l0;j++)
						{
							if(powerfileflag == 0)
		   	   			apowers = powerisotherm(akt);																	// apower is the initial power spectrum
         	  			else if(powerfileflag == 1)
								apowers = powerI;

         	  			apowers=apowers/akt;
           			
           				ckj=apowers*dliso[j]*dliso[j]*dak1[k];
           				cliso[j][in]=cliso[j][in]+ckj;
           			
        					cpkj=apowers*dpliso[j]*dpliso[j]*dak1[k];															// Calculating Cl's 	
           				cpliso[j][in]=cpliso[j][in]+cpkj;																	// Zaldarriaga PhD Thesis 3.25	
            		
           				cckj=apowers*dliso[j]*dpliso[j]*dak1[k];
        				 	ccliso[j][in]=ccliso[j][in]+cckj;
  		     				}

		     			for(int j=1;j<=l0;j++)
						{
							if(powerfileflag == 0)
		   	   			apowers = powercross(akt);																	// apower is the initial power spectrum
         	  			else if(powerfileflag == 1)
								apowers = powerC;

         	  			apowers=apowers/akt;
           			
           				ckj=apowers*dladb[j]*dliso[j]*dak1[k];
           				clcross[j][in]=clcross[j][in]+ckj;
           			
        					cpkj=apowers*dpladb[j]*dpliso[j]*dak1[k];															// Calculating Cl's 	
           				cplcross[j][in]=cplcross[j][in]+cpkj;																	// Zaldarriaga PhD Thesis 3.25	
            		
           				cckj=apowers*dladb[j]*dpliso[j]*dak1[k];
        				 	cclcross[j][in]=cclcross[j][in]+cckj;
							//printf("%e \t%e \t%e\n",cladb[j][in],cliso[j][in],clcross[j][in]);
		     				}
						}
					}
				}				// IF OF ITFLAG != 2 IE SCALER PART LOOP ENDS HERE
			/***************************************************************************************************************************/
			/*     Begin l and  time-loop to integrate tensor perturbations.																				*/
			/*     Finding the ranges of integration.																												*/
			/***************************************************************************************************************************/
			
		   if (itflag!=0)
   		{
   		   for(int j=1;j<=l0;j++)
   		   {
	     			xlim=0.05*l[j];
        			xlim=max(xlim,xlimmin);
        			xlim=l[j]-xlim;
        			xlmax1=80.0*l[j];
        			tmin=tau0-xlmax1/akt;
        			tmax=tau0-xlim/akt;
        			tmax=min(tau0,tmax);
        			tmin=max(atau0[2],tmin);
				
     			   if (tmax < atau0[2]) goto MOVETOB;
				
     			   if (zri==0.0)
     			   {
	        			if (tmin < taurend)
	        				nstart1=2;
        				 else
           				nstart1=n1+int(log(tmin/taurend)/dlntau0);
				
   		         if (tmax<taurend)
        			   {
	        			   nstop1=n1;
           				nstop1a=n1;
           				}
           			else
           			{
              			nstop1=n1+(int)(log(tmax/taurend)/dlntau0);
           				nstop1=minint(nstop1,nstpt);
           	
           				if ((akt*dtau1[nstop1])>1)
           				{
	           				nstop1a=maxint(n1,n1-(int)(log(akt*dlntau0*taurend)/dlntau0));
           					}
           				else
           				{
           					nstop1a=nstop1;
           					}
        					}
        				}	
     				
     				else
     				{
	     				if (tmin < taurend)
        				{
	        				nstart1=2;
           				nstart2=j2ri1;
        	 				}
         	
        				else
        				{
           				if (tmin < atau0[nriend])
           				{   
              				nstart1=minint(n1+int(log(tmin/taurend)/dlntau0),j2ri1);
           					nstart2=j2ri1;
           					}
           				else
           				{
           				 	nstart1=nstep+1;
   	     					nstart2=nriend+int(log(tmin/atau0[nriend]/dlntau0));
        						}
        					}

	     			   if (tmax < taurend)
   	  			   {
         			   nstop1=n1;
        	   			nstop2=0;
        					}
        	
        				else
        				{
        					if (tmax < atau0[j2ri1])
        	 				{
              				nstop1=n1+(int)(log(tmax/taurend)/dlntau0);
           					nstop2=0;
           					}
           	
        				 	else
           				{
              				nstop1=j2ri1-1;
              				nstop2=maxint(nriend,nriend+(int)(log(tmax/atau0[nriend])/dlntau0));

               			nstop2=minint(nstop2,nstpt);
	               
		               	if ((akt*dtau1[nstop2]) > 1)
        			      	{
              					nstop2a=maxint(nriend,nriend-(int)(log(akt*dlntau0*atau0[nriend])/dlntau0));
              					}
              
           			   	else
              				{
              					nstop2a=nstop2;
           						}
            				}
					
           				nstop1=minint(nstop1,j2ri1-1);
          		
           				if ((akt*dtau1[nstop1]) > 1)
           					nstop1a=maxint(n1,n1-int(log(akt*dlntau0*taurend)/dlntau0));
           		
           				else
           					nstop1a=nstop1;
        					}
						}
					/***************************************************************************************************************************/
					/*    Integration before reionization																													*/
					/*		Interpolating jls at points where the sources are recorded.																				*/
					/***************************************************************************************************************************/			
         
              	for(int i=nstart1; i<=nstop1a; i++)
              	{
        	     	   xf=akt*(tau0-atau0[i]);
        	      
                 	if ((j==1) && (xf<1.0) && (xf>1.0e-7))
                 		ajl0=(3.0*(sin(xf)/xf-cos(xf))/xf - sin(xf))/xf;
                 	else
                 	{
                 		m2=mxx[i];
                    	h2=xx[m2+1]-xx[m2];
                    	a2=(xx[m2+1]-xf)/h2;
                    	b2=(xf-xx[m2])/h2;
                    	ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
                    	}
                             
                 	dtl2[j]=dtl2[j]+st2[i]*ajl0*dtau2[i];
                 	dtel2[j]=dtel2[j]+ste2[i]*ajl0*dtau2[i];
                 	dtbl2[j]=dtbl2[j]+stb2[i]*ajl0*dtau2[i];
              		}
         		
        			/***************************************************************************************************************************/
        			/*     Tensor case.																																			*/
        			/*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
        			/***************************************************************************************************************************/		
         		
        	      for(int i=nstop1a+1; i<=nstop1; i++)
              	{
        	      	xf=akt*(tau0-atau0[i]);
                 	m2=mxx[i];
                 	xi=xf-akt*dtau1[i];
                 	m1=mxx[i+1];
                 	dtdt1=st2[i];
                 	dtdt2=st2[i+1];
                 	dtde1=ste2[i];
                 	dtde2=ste2[i+1];
                	dtdb1=stb2[i];
                 	dtdb2=stb2[i+1];
                  
                 	for(int lx=m1+2;lx<=m2+1;lx++)
                 	{
        	         	x=xx[lx];
                    	dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2;
                    	dtde=(dtde1-dtde2)*(x-xi)/(xf-xi)+dtde2;
                    	dtdb=(dtdb1-dtdb2)*(x-xi)/(xf-xi)+dtdb2;
                    	dtl3[j]=dtl3[j]+ajl[lx][j]*dtdt*dxx[lx];
                    	dtel3[j]=dtel3[j]+ajl[lx][j]*dtde*dxx[lx];
                    	dtbl3[j]=dtbl3[j]+ajl[lx][j]*dtdb*dxx[lx];
                    	}
                 	}
                  
        			/***************************************************************************************************************************/
        			/*    Integration after reionization																													*/
        			/*		Interpolating jls at points where the sources are recorded.																				*/
        			/***************************************************************************************************************************/
         
        			if (zri!=0.0) 
        			{
        
                 	for(int i=nstart2;i <= nstop2a;i++)
                 	{
        	        		xf=akt*(tau0-atau0[i]);
        
                    	if ((j==1) && (xf<1.0) && (xf>1.0e-7))
                    		ajl0=(3.0*(sin(xf)/xf-cos(xf))/xf- sin(xf))/xf;
                    	else
                    	{
        	            	m2=mxx[i];
                       	h2=xx[m2+1]-xx[m2];
                       	a2=(xx[m2+1]-xf)/h2;
                       	b2=(xf-xx[m2])/h2;
                       	ajl0=a2*ajl[m2][j]+b2*ajl[m2+1][j]+((a2*a2*a2-a2)*ajlpr[m2][j]+(b2*b2*b2-b2)*ajlpr[m2+1][j])*(h2*h2)/6.0;
                       	}
                       
                    	dtl2[j]=dtl2[j]+st2[i]*ajl0*dtau2[i];
                    	dtel2[j]=dtel2[j]+ste2[i]*ajl0*dtau2[i];
                    	dtbl2[j]=dtbl2[j]+stb2[i]*ajl0*dtau2[i];
                    	}
        
        
        				/***************************************************************************************************************************/
        				/*     Tensor case.																																			*/
        				/*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
        				/***************************************************************************************************************************/		
        				
        				for(int i=nstop2a+1;i<=nstop2;i++)
        				{
        					xi=xf-akt*dtau1[i];
                    	xf=akt*(tau0-atau0[i]);
                    	m2=mxx[i];
                    	m1=mxx[i+1];
                    	dtdt1=st2[i];
                    	dtdt2=st2[i+1];
                    	dtde1=ste2[i];
                    	dtde2=ste2[i+1];
                    	dtdb1=stb2[i];
                    	dtdb2=stb2[i+1];
                    
                    	for(int lx=m1+2;lx <= m2+1;lx++)
                    	{
        	            	x=xx[lx];
                       	dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2;
                       	dtde=(dtde1-dtde2)*(x-xi)/(xf-xi)+dtde2;
                       	dtdb=(dtdb1-dtdb2)*(x-xi)/(xf-xi)+dtdb2;
                       	dtl3[j]=dtl3[j]+ajl[lx][j]*dtdt*dxx[lx];
                       	dtel3[j]=dtel3[j]+ajl[lx][j]*dtde*dxx[lx];
                       	dtbl3[j]=dtbl3[j]+ajl[lx][j]*dtdb*dxx[lx];
                       	}
                    	}
                 	}
					                  
               dtl[j]=dtl2[j]+dtl3[j]/akt;
               dtel[j]=dtel2[j]+dtel3[j]/akt;
               dtbl[j]=dtbl2[j]+dtbl3[j]/akt;
					}

MOVETOB:			         
				/***************************************************************************************************************************/
				/*     Tensor case.																																			*/
				/*		 Finally calculate Cl by intigrating over K's . ps: l*(l-1) will be done in the next step										*/
				/***************************************************************************************************************************/
					
         	for(int in=1;in<=nn;in++)																				
         	{
         		for(int j=1;j<=l0;j++)
         		{
         			powertflat(akt,in,&apowert);																	// apower is the initial power spectrum
              		apowert=apowert/akt;	
         
         			ctkj=apowert*dtl[j]*dtl[j]*dak1[k];															// Calculating Cl's 		
               	ctl[j][in]=ctl[j][in]+ctkj;																	// Zaldarriaga PhD Thesis 3.25	
         	
               	ctekj=apowert*dtel[j]*dtel[j]*dak1[k];
           	   	ctel[j][in]=ctel[j][in]+ctekj;
         
              		ctbkj=apowert*dtbl[j]*dtbl[j]*dak1[k];
              		ctbl[j][in]=ctbl[j][in]+ctbkj;
         
              		ctckj=apowert*dtl[j]*dtel[j]*dak1[k];
                	ctcl[j][in]=ctcl[j][in]+ctckj;
         			}
         		}	
         	}  		// IF ITFLAG !=0 ENDS HERE    
         }				// K LOOP ENDS HERE   
			
      // (P18) : Calculate the Cl's and then calculate the primes at all these points  
         
		/******************************************************************************************************************************/
		/*  Final calculations for CMB output.																														*/
		/******************************************************************************************************************************/
					
      for(int j=1;j<=l0;j++)
      	xl[j]=l[j];
		
		/******************************************************************************************************************************/
	   /*    Scalar case																																					*/	
      /******************************************************************************************************************************/
      if (itflag!=2)
      {
			if(twofield == 0)
			{
				for(int in=1; in<=nn;in++)
   	   	{
   	     	   for(int j=1;j<=l0;j++)
   	        	{
						ctnorm=l[j]*(l[j]-1.0)*(l[j]+1)*(l[j]+2);
   	            cl[j][in]=2.0*cl[j][in]*l[j]*(l[j]+1);
   	            cpl[j][in]=2.0*ctnorm*cpl[j][in]*l[j]*(l[j]+1);
   	            ccl[j][in]=2.0*sqrt(ctnorm)*ccl[j][in]*l[j]*(l[j]+1);
   	            ckkl[j][in]=2.0*ckkl[j][in]*l[j]*(l[j]+1);
   		         ctkl[j][in]=2.0*ctkl[j][in]*l[j]*(l[j]+1);
						//printf("%e\n",cl[j][in]);
     					}
					}
				}
			
			if(twofield == 1)
			{
				for(int in=1; in<=nn;in++)
   	   	{
   	     	   for(int j=1;j<=l0;j++)
   	        	{
						ctnorm=l[j]*(l[j]-1.0)*(l[j]+1)*(l[j]+2);
   	            cladb[j][in]=2.0*cladb[j][in]*l[j]*(l[j]+1);
   	            cpladb[j][in]=2.0*ctnorm*cpladb[j][in]*l[j]*(l[j]+1);
   	            ccladb[j][in]=2.0*sqrt(ctnorm)*ccladb[j][in]*l[j]*(l[j]+1);
     					}
					}

				for(int in=1; in<=nn;in++)
   	   	{
   	     	   for(int j=1;j<=l0;j++)
   	        	{
						ctnorm=l[j]*(l[j]-1.0)*(l[j]+1)*(l[j]+2);
   	            cliso[j][in]=2.0*cliso[j][in]*l[j]*(l[j]+1);
   	            cpliso[j][in]=2.0*ctnorm*cpliso[j][in]*l[j]*(l[j]+1);
   	            ccliso[j][in]=2.0*sqrt(ctnorm)*ccliso[j][in]*l[j]*(l[j]+1);
     					}
					}
				
				for(int in=1; in<=nn;in++)
   	   	{
   	     	   for(int j=1;j<=l0;j++)
   	        	{
						ctnorm=l[j]*(l[j]-1.0)*(l[j]+1)*(l[j]+2);
   	            clcross[j][in]=2.0*clcross[j][in]*l[j]*(l[j]+1);
   	            cplcross[j][in]=2.0*ctnorm*cplcross[j][in]*l[j]*(l[j]+1);
   	            cclcross[j][in]=2.0*sqrt(ctnorm)*cclcross[j][in]*l[j]*(l[j]+1);
     					}
					}
				
				for(int in=1; in<=nn;in++)
   	   	{
   	     	   for(int j=1;j<=l0;j++)
   	        	{
   	            //cltotal[j][in] =-(4.0/15)*clcross[j][in]  + (1.0/9)*cladb[j][in]  + (4.0/25)*cliso[j][in];
   	            //cpltotal[j][in]=-(4.0/15)*cplcross[j][in] + (1.0/9)*cpladb[j][in] + (4.0/25)*cpliso[j][in];
   	            //ccltotal[j][in]=-(4.0/15)*cclcross[j][in] + (1.0/9)*ccladb[j][in] + (4.0/25)*ccliso[j][in];
 
   	            cltotal[j][in] =2*clcross[j][in]  + cladb[j][in]  + cliso[j][in];
   	            cpltotal[j][in]=2*cplcross[j][in] + cpladb[j][in] + cpliso[j][in];
   	            ccltotal[j][in]=2*cclcross[j][in] + ccladb[j][in] + ccliso[j][in];
    					
						}
					}
				}

			/***************************************************************************************************************************/
		   /*   Making the interpolation tables to get other l-values.																						*/	
	      /***************************************************************************************************************************/               

         llo=1;
         cllo=1.0e30;                               
         clhi=1.0e30;
            
			if(twofield == 0)
			{ 	        
				for(int in=1;in<=nn;in++)
	         {
					for(int jj=1;jj<=lmax+1;jj++)
	 				{
		 				bpcl[jj]=cl[jj][in];
		 				bpcpl[jj]=cpl[jj][in];
		 				bpccl[jj]=ccl[jj][in];
		 				bpckkl[jj]=ckkl[jj][in];
		 				bpctkl[jj]=ctkl[jj][in];
		 				}
		         
		         spline(xl,bpcl,l0,cllo,clhi,bpclpr);
	            spline(xl,bpcpl,l0,cllo,clhi,bpcplpr);
	            spline(xl,bpccl,l0,cllo,clhi,bpcclpr);
	            spline(xl,bpckkl,l0,cllo,clhi,bpckklpr);
	            spline(xl,bpctkl,l0,cllo,clhi,bpctklpr);
		
  					for(int jj=1;jj<=lmax+1;jj++)
 					{
		 				clpr[jj][in]=bpclpr[jj];
		 				cplpr[jj][in]=bpcplpr[jj];
		 				cclpr[jj][in]=bpcclpr[jj];
		 				ckklpr[jj][in]=bpckklpr[jj];
		 				ctklpr[jj][in]=bpctklpr[jj];
		 				}
      	      }
				}


			if(twofield == 1)
			{ 	        
				for(int in=1;in<=nn;in++)
	         {
					for(int jj=1;jj<=lmax+1;jj++)
	 				{
		 				bpcl[jj]=cladb[jj][in];
		 				bpcpl[jj]=cpladb[jj][in];
		 				bpccl[jj]=ccladb[jj][in];
		 				}
		         
		         spline(xl,bpcl,l0,cllo,clhi,bpclpr);
	            spline(xl,bpcpl,l0,cllo,clhi,bpcplpr);
	            spline(xl,bpccl,l0,cllo,clhi,bpcclpr);
		
  					for(int jj=1;jj<=lmax+1;jj++)
 					{
		 				clpradb[jj][in]=bpclpr[jj];
		 				cplpradb[jj][in]=bpcplpr[jj];
		 				cclpradb[jj][in]=bpcclpr[jj];
		 				}
      	      }

				for(int in=1;in<=nn;in++)
	         {
					for(int jj=1;jj<=lmax+1;jj++)
	 				{
		 				bpcl[jj]=cliso[jj][in];
		 				bpcpl[jj]=cpliso[jj][in];
		 				bpccl[jj]=ccliso[jj][in];
		 				}
		         
		         spline(xl,bpcl,l0,cllo,clhi,bpclpr);
	            spline(xl,bpcpl,l0,cllo,clhi,bpcplpr);
	            spline(xl,bpccl,l0,cllo,clhi,bpcclpr);
		
  					for(int jj=1;jj<=lmax+1;jj++)
 					{
		 				clpriso[jj][in]=bpclpr[jj];
		 				cplpriso[jj][in]=bpcplpr[jj];
		 				cclpriso[jj][in]=bpcclpr[jj];
		 				}
      	      }

				for(int in=1;in<=nn;in++)
	         {
					for(int jj=1;jj<=lmax+1;jj++)
	 				{
		 				bpcl[jj]=clcross[jj][in];
		 				bpcpl[jj]=cplcross[jj][in];
		 				bpccl[jj]=cclcross[jj][in];
		 				}
		         
		         spline(xl,bpcl,l0,cllo,clhi,bpclpr);
	            spline(xl,bpcpl,l0,cllo,clhi,bpcplpr);
	            spline(xl,bpccl,l0,cllo,clhi,bpcclpr);
		
  					for(int jj=1;jj<=lmax+1;jj++)
 					{
		 				clprcross[jj][in]=bpclpr[jj];
		 				cplprcross[jj][in]=bpcplpr[jj];
		 				cclprcross[jj][in]=bpcclpr[jj];
		 				}
      	      }

				for(int in=1;in<=nn;in++)
	         {
					for(int jj=1;jj<=lmax+1;jj++)
	 				{
		 				bpcl[jj]=cltotal[jj][in];
		 				bpcpl[jj]=cpltotal[jj][in];
		 				bpccl[jj]=ccltotal[jj][in];
		 				}
		         
		         spline(xl,bpcl,l0,cllo,clhi,bpclpr);
	            spline(xl,bpcpl,l0,cllo,clhi,bpcplpr);
	            spline(xl,bpccl,l0,cllo,clhi,bpcclpr);
		
  					for(int jj=1;jj<=lmax+1;jj++)
 					{
		 				clprtotal[jj][in]=bpclpr[jj];
		 				cplprtotal[jj][in]=bpcplpr[jj];
		 				cclprtotal[jj][in]=bpcclpr[jj];
		 				}
      	      }

				}
	      }
		
		/***************************************************************************************************************************/
		/*     Tensor Case																																			*/
		/***************************************************************************************************************************/
		if (itflag!=0)
      {
			/************************************************************************************************************************/	      
			/*    Normalization																																		*/
			/************************************************************************************************************************/
			for(int j=1;j<=l0;j++)
   	   {
	   	   ctnorm=l[j]*(l[j]-1.0)*(l[j]+1)*(l[j]+2);
         	for(int in=1;in<=nn;in++)
            {
	            ctl[j][in]=ctnorm*ctl[j][in]*(l[j]*(l[j]+1))/8.0;
               ctel[j][in]=ctel[j][in]*(l[j]*(l[j]+1))/8.0;
               ctbl[j][in]=ctbl[j][in]*(l[j]*(l[j]+1))/8.0;
               ctcl[j][in]=sqrt(ctnorm)*ctcl[j][in]*l[j]*(l[j]+1)/8.0;
 					}
 				}

			/************************************************************************************************************************/	      
			/*		Making the interpolation tables to get other l-values.																				*/
			/************************************************************************************************************************/	      

			cllo=1.0e30;
         clhi=1.0e30;
			               
         for(int in=1;in<=nn;in++)
         {
	         for(int jj=1;jj<=lmax+1;jj++)
 				{
	 				bpctl[jj]=ctl[jj][in];
	 				bpctel[jj]=ctel[jj][in];
	 				bpctbl[jj]=ctbl[jj][in];
	 				bpctcl[jj]=ctcl[jj][in];
	 				}

            spline(xl,bpctl,l0,cllo,clhi,bpctlpr);
            spline(xl,bpctel,l0,cllo,clhi,bpctelpr);
            spline(xl,bpctbl,l0,cllo,clhi,bpctblpr);
            spline(xl,bpctcl,l0,cllo,clhi,bpctclpr);

 				for(int jj=1;jj<=lmax+1;jj++)
		 		{
			 		ctlpr[jj][in]=bpctlpr[jj];
			 		ctelpr[jj][in]=bpctelpr[jj];
	 				ctblpr[jj][in]=bpctblpr[jj];
	 				ctclpr[jj][in]=bpctclpr[jj];
	 				}
    			}
		   }				// IF ITFLAG!=1 ENDS HERE

		/************************************************************************************************************************/	      
		/*     Calculating Cls for every l.																													*/
		/************************************************************************************************************************/	      

		// (P20) : Interpolate the Cls for all the l's

      for(int in=1;in<=nn;in++)
      {
	      llo=1;
         for(int il=2;il<=l[l0];il++)
         {
	         xi=il;
                  
            if ((xi > xl[llo+1]) && (llo < l0))
            	llo=llo+1;
                                 
            lhi=llo+1;
            ho=xl[lhi]-xl[llo];
            a0=(xl[lhi]-xi)/ho;
            b0=(xi-xl[llo])/ho;
                  
            if (itflag!=2)
            {
					if(twofield == 0)
					{
		            clint=a0*cl[llo][in]+b0*cl[lhi][in]+((pow(a0,3)-a0)*clpr[llo][in]+(pow(b0,3)-b0)*clpr[lhi][in])*ho*ho /6.0;
   	            cplint=a0*cpl[llo][in]+b0*cpl[lhi][in]+((pow(a0,3)-a0)*cplpr[llo][in]+(pow(b0,3)-b0)*cplpr[lhi][in])*ho*ho /6.0;
      	         cclint=a0*ccl[llo][in]+b0*ccl[lhi][in]+((pow(a0,3)-a0)*cclpr[llo][in]+(pow(b0,3)-b0)*cclpr[lhi][in])*ho*ho /6.0;
         	      ckklint=a0*ckkl[llo][in]+b0*ckkl[lhi][in]+((pow(a0,3)-a0)*ckklpr[llo][in]+(pow(b0,3)-b0)*ckklpr[lhi][in])*ho*ho /6.0;
            	   ctklint=a0*ctkl[llo][in]+b0*ctkl[lhi][in]+((pow(a0,3)-a0)*ctklpr[llo][in]+(pow(b0,3)-b0)*ctklpr[lhi][in])*ho*ho /6.0;
               	}
					if(twofield == 1)
					{
		            clintadb=a0*cladb[llo][in]+b0*cladb[lhi][in]+((pow(a0,3)-a0)*clpradb[llo][in]+(pow(b0,3)-b0)*clpradb[lhi][in])*ho*ho /6.0;
   	            cplintadb=a0*cpladb[llo][in]+b0*cpladb[lhi][in]+((pow(a0,3)-a0)*cplpradb[llo][in]+(pow(b0,3)-b0)*cplpradb[lhi][in])*ho*ho /6.0;
      	         cclintadb=a0*ccladb[llo][in]+b0*ccladb[lhi][in]+((pow(a0,3)-a0)*cclpradb[llo][in]+(pow(b0,3)-b0)*cclpradb[lhi][in])*ho*ho /6.0;

		            clintiso=a0*cliso[llo][in]+b0*cliso[lhi][in]+((pow(a0,3)-a0)*clpriso[llo][in]+(pow(b0,3)-b0)*clpriso[lhi][in])*ho*ho /6.0;
   	            cplintiso=a0*cpliso[llo][in]+b0*cpliso[lhi][in]+((pow(a0,3)-a0)*cplpriso[llo][in]+(pow(b0,3)-b0)*cplpriso[lhi][in])*ho*ho /6.0;
      	         cclintiso=a0*ccliso[llo][in]+b0*ccliso[lhi][in]+((pow(a0,3)-a0)*cclpriso[llo][in]+(pow(b0,3)-b0)*cclpriso[lhi][in])*ho*ho /6.0;

		            clintcross=a0*clcross[llo][in]+b0*clcross[lhi][in]+((pow(a0,3)-a0)*clprcross[llo][in]+(pow(b0,3)-b0)*clprcross[lhi][in])*ho*ho /6.0;
   	            cplintcross=a0*cplcross[llo][in]+b0*cplcross[lhi][in]+((pow(a0,3)-a0)*cplprcross[llo][in]+(pow(b0,3)-b0)*cplprcross[lhi][in])*ho*ho /6.0;
      	         cclintcross=a0*cclcross[llo][in]+b0*cclcross[lhi][in]+((pow(a0,3)-a0)*cclprcross[llo][in]+(pow(b0,3)-b0)*cclprcross[lhi][in])*ho*ho /6.0;

		            clinttotal=a0*cltotal[llo][in]+b0*cltotal[lhi][in]+((pow(a0,3)-a0)*clprtotal[llo][in]+(pow(b0,3)-b0)*clprtotal[lhi][in])*ho*ho /6.0;
   	            cplinttotal=a0*cpltotal[llo][in]+b0*cpltotal[lhi][in]+((pow(a0,3)-a0)*cplprtotal[llo][in]+(pow(b0,3)-b0)*cplprtotal[lhi][in])*ho*ho /6.0;
      	         cclinttotal=a0*ccltotal[llo][in]+b0*ccltotal[lhi][in]+((pow(a0,3)-a0)*cclprtotal[llo][in]+(pow(b0,3)-b0)*cclprtotal[lhi][in])*ho*ho /6.0;
						}
					} 
                  
            if (itflag!=0)
            {
	            ctlint=a0*ctl[llo][in]+b0*ctl[lhi][in]+((pow(a0,3)-a0)*ctlpr[llo][in]+(pow(b0,3)-b0)*ctlpr[lhi][in])*ho*ho/6.0;
               ctelint=a0*ctel[llo][in]+b0*ctel[lhi][in]+((pow(a0,3)-a0)*ctelpr[llo][in]+(pow(b0,3)-b0)*ctelpr[lhi][in])*ho*ho/6.0;
               ctblint=a0*ctbl[llo][in]+b0*ctbl[lhi][in]+((pow(a0,3)-a0)*ctblpr[llo][in]+(pow(b0,3)-b0)*ctblpr[lhi][in])*ho*ho/6.0;
               ctclint=a0*ctcl[llo][in]+b0*ctcl[lhi][in]+((pow(a0,3)-a0)*ctclpr[llo][in]+(pow(b0,3)-b0)*ctclpr[lhi][in])*ho*ho/6.0;
               } 

   
           	if (itflag!=2)
           	{
					if(twofield == 0)
					{
		         	clts[il][in]=clint/cl[1][in];
   	           	cles[il][in]=cplint/cl[1][in];
      	        	clcs[il][in]=cclint/cl[1][in];
         	     	clkk[il][in]=ckklint/cl[1][in];
            	  	cltk[il][in]=ctklint/cl[1][in];

               	}
					if(twofield == 1)
					{
						/*
		         	cltsadb[il][in]=clintadb/cladb[1][in];
   	           	clesadb[il][in]=cplintadb/cladb[1][in];
      	        	clcsadb[il][in]=cclintadb/cladb[1][in];

		         	cltsiso[il][in]=clintiso/cliso[1][in];
   	           	clesiso[il][in]=cplintiso/cliso[1][in];
      	        	clcsiso[il][in]=cclintiso/cliso[1][in];

		         	cltscross[il][in]=clintcross/clcross[1][in];
   	           	clescross[il][in]=cplintcross/clcross[1][in];
      	        	clcscross[il][in]=cclintcross/clcross[1][in];
						*/

		         	cltsadb[il][in]=clintadb/cltotal[1][in];
   	           	clesadb[il][in]=cplintadb/cltotal[1][in];
      	        	clcsadb[il][in]=cclintadb/cltotal[1][in];

		         	cltsiso[il][in]=clintiso/cltotal[1][in];
   	           	clesiso[il][in]=cplintiso/cltotal[1][in];
      	        	clcsiso[il][in]=cclintiso/cltotal[1][in];

		         	cltscross[il][in]=clintcross/cltotal[1][in];
   	           	clescross[il][in]=cplintcross/cltotal[1][in];
      	        	clcscross[il][in]=cclintcross/cltotal[1][in];

		         	cltstotal[il][in]=clinttotal/cltotal[1][in];
   	           	clestotal[il][in]=cplinttotal/cltotal[1][in];
      	        	clcstotal[il][in]=cclinttotal/cltotal[1][in];
						}
					}

           	if (itflag!=0)
           	{
	          	cltt[il][in]=ctlint/ctl[1][in];
              	clet[il][in]=ctelint/ctl[1][in];
              	clbt[il][in]=ctblint/ctl[1][in];
              	clct[il][in]=ctclint/ctl[1][in];
              	}
   			}
			   					
         if (itflag!=2)
			{
				if(twofield == 0)
	         	clts[1][in]=cl[1][in];
				else if(twofield == 1)
				{
	         	cltsadb[1][in] = cltotal[1][in];
	         	cltsiso[1][in] = cltotal[1][in];
	         	cltscross[1][in] = cltotal[1][in];
	         	cltstotal[1][in] = cltotal[1][in];					
					}
				}
               	
         if (itflag!=0)
         	cltt[1][in]=ctl[1][in];
         }
		}	
	}	







void COBEnormalizetwofield()
{

   const double fourpi=4.0*3.14159265;
	double xlog10,h,xlnh,hc,curv,r;
	double c10,d1,d2,d3,d4,d5,d6,d7,x1,x2,x3,x4,x5,x6,x7;
	double sy,s,sx,sxy,sxx,delt,d1pr,d1ppr;
	FILE *fp;
        fp=fopen("cltwofieldx.dat","w");
	double d2norm;								// Needed for lensing calculation	
		
   xlog10=log(10.0);
   h=h0/100.0;
   xlnh=log(h);

	// Curvature radius
	if(fabs(omegak)>1.0e-3)
  	{
		hc=2.998e5/h0;                      // h*c. c = 2.338*e5 km/s. h = km/sec/pc 
      curv=-omegak/(hc*hc);					// Gaussian curvature	
      r=1.0/sqrt(fabs(curv));				   // Radious of curvature	
      }

	// COBE normalization
	// fit the spectrum to a quadratic around C_10 with equal weights in logl
	/*
	for(int in=1;in<=nn;in++)
	{
		c10=cltsadb[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,cltsadb[10][in]);
      d1=(cltsadb[l[2]][in])/c10-1.0;
      d2=(cltsadb[l[3]][in])/c10-1.0;
      d3=(cltsadb[l[5]][in])/c10-1.0;
      d4=(cltsadb[l[7]][in])/c10-1.0;
      d5=(cltsadb[l[10]][in])/c10-1.0;
      d6=(cltsadb[l[11]][in])/c10-1.0;
      d7=(cltsadb[l[12]][in])/c10-1.0;
		
      x1=log(1.0*l[2])/xlog10-1.0;
      x2=log(1.0*l[3])/xlog10-1.0;
      x3=log(1.0*l[5])/xlog10-1.0;
      x4=log(1.0*l[7])/xlog10-1.0;
      x5=log(1.0*l[10])/xlog10-1.0;
      x6=log(1.0*l[11])/xlog10-1.0;
      x7=log(1.0*l[12])/xlog10-1.0;

		// Quadratic List Square fit      
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);
      
      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;
		
		printf("%e %e %e\n",delt,d1pr,d1ppr);
		
		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;
		
		// C_l normalization; output l(l+1)C_l/twopi
		
		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);		
		for(int il=2;il<=l[l0];il++)
		{
			printf("clts = %e \n",cltsadb[il][in]);
			cltsadb[il][in]=cltsadb[il][in]*c10;
	      clesadb[il][in]=clesadb[il][in]*c10;
   	   clcsadb[il][in]=clcsadb[il][in]*c10;
   	   }
		}
	*/
	/*
	for(int in=1;in<=nn;in++)
	{
		c10=cltsiso[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(cltsiso[l[2]][in])/c10-1.0;
      d2=(cltsiso[l[3]][in])/c10-1.0;
      d3=(cltsiso[l[5]][in])/c10-1.0;
      d4=(cltsiso[l[7]][in])/c10-1.0;
      d5=(cltsiso[l[10]][in])/c10-1.0;
      d6=(cltsiso[l[11]][in])/c10-1.0;
      d7=(cltsiso[l[12]][in])/c10-1.0;
		
      x1=log(1.0*l[2])/xlog10-1.0;
      x2=log(1.0*l[3])/xlog10-1.0;
      x3=log(1.0*l[5])/xlog10-1.0;
      x4=log(1.0*l[7])/xlog10-1.0;
      x5=log(1.0*l[10])/xlog10-1.0;
      x6=log(1.0*l[11])/xlog10-1.0;
      x7=log(1.0*l[12])/xlog10-1.0;

		// Quadratic List Square fit      
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);
      
      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;
		
		printf("%e %e %e\n",delt,d1pr,d1ppr);
		
		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;
		
		// C_l normalization; output l(l+1)C_l/twopi
		
		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);		
		for(int il=2;il<=l[l0];il++)
		{
			printf("clts = %e \n",cltsiso[il][in]);
			cltsiso[il][in]=cltsiso[il][in]*c10;
	      clesiso[il][in]=clesiso[il][in]*c10;
   	   clcsiso[il][in]=clcsiso[il][in]*c10;
   	   }
		}
	*/
	/*
	for(int in=1;in<=nn;in++)
	{
		c10=cltscross[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(cltscross[l[2]][in])/c10-1.0;
      d2=(cltscross[l[3]][in])/c10-1.0;
      d3=(cltscross[l[5]][in])/c10-1.0;
      d4=(cltscross[l[7]][in])/c10-1.0;
      d5=(cltscross[l[10]][in])/c10-1.0;
      d6=(cltscross[l[11]][in])/c10-1.0;
      d7=(cltscross[l[12]][in])/c10-1.0;
		
      x1=log(1.0*l[2])/xlog10-1.0;
      x2=log(1.0*l[3])/xlog10-1.0;
      x3=log(1.0*l[5])/xlog10-1.0;
      x4=log(1.0*l[7])/xlog10-1.0;
      x5=log(1.0*l[10])/xlog10-1.0;
      x6=log(1.0*l[11])/xlog10-1.0;
      x7=log(1.0*l[12])/xlog10-1.0;

		// Quadratic List Square fit      
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);
      
      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;
		
		printf("%e %e %e\n",delt,d1pr,d1ppr);
		
		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;
		
		// C_l normalization; output l(l+1)C_l/twopi
		
		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);		
		for(int il=2;il<=l[l0];il++)
		{
			printf("clts = %e \n",clts[il][in]);
			cltscross[il][in]=cltscross[il][in]*c10;
	      clescross[il][in]=clescross[il][in]*c10;
   	   clcscross[il][in]=clcscross[il][in]*c10;
   	   }
		}
	*/
	for(int in=1;in<=nn;in++)
	{
		c10=cltstotal[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(cltstotal[l[2]][in])/c10-1.0;
      d2=(cltstotal[l[3]][in])/c10-1.0;
      d3=(cltstotal[l[5]][in])/c10-1.0;
      d4=(cltstotal[l[7]][in])/c10-1.0;
      d5=(cltstotal[l[10]][in])/c10-1.0;
      d6=(cltstotal[l[11]][in])/c10-1.0;
      d7=(cltstotal[l[12]][in])/c10-1.0;
		
      x1=log(1.0*l[2])/xlog10-1.0;
      x2=log(1.0*l[3])/xlog10-1.0;
      x3=log(1.0*l[5])/xlog10-1.0;
      x4=log(1.0*l[7])/xlog10-1.0;
      x5=log(1.0*l[10])/xlog10-1.0;
      x6=log(1.0*l[11])/xlog10-1.0;
      x7=log(1.0*l[12])/xlog10-1.0;

		// Quadratic List Square fit      
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);
      
      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;
		
		printf("%e %e %e\n",delt,d1pr,d1ppr);
		
		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;
		
		// C_l normalization; output l(l+1)C_l/twopi
		
		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);		
		for(int il=2;il<=l[l0];il++)
		{
			printf("clts = %e \n",cltstotal[il][in]);

			cltsadb[il][in]=(1/9.0)*cltsadb[il][in]*c10;
	      clesadb[il][in]=(1/9.0)*clesadb[il][in]*c10;
   	   clcsadb[il][in]=(1/9.0)*clcsadb[il][in]*c10;

			cltsiso[il][in]=(4.0/25)*cltsiso[il][in]*c10;
	      clesiso[il][in]=(4.0/25)*clesiso[il][in]*c10;
   	   clcsiso[il][in]=(4.0/25)*clcsiso[il][in]*c10;

			cltscross[il][in]=-(4.0/15)*cltscross[il][in]*c10;
	      clescross[il][in]=-(4.0/15)*clescross[il][in]*c10;
   	   clcscross[il][in]=-(4.0/15)*clcscross[il][in]*c10;

			cltstotal[il][in]=cltstotal[il][in]*c10;
	      clestotal[il][in]=clestotal[il][in]*c10;
   	   clcstotal[il][in]=clcstotal[il][in]*c10;
           
           fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e \n",cltsadb[il][in],clesadb[il][in],clcsadb[il][in],cltsiso[il][in],clesiso[il][in],clcsiso[il][in],cltscross[il][in],clescross[il][in],clcscross[il][in],cltstotal[il][in],clestotal[il][in],clcstotal[il][in]);
   	   }
		}
	}




/*********************************************************************************/
/* Normalization is requiredfor matching the data with the observed data set     */
/* Otherwise different theory will give different amplitude of the spectrum      */
/* and it would be difficult to find out the accuricy of the theory              */
/*                                                                               */
/* This part is based on the paper: 													      */
/*      The Four-Year COBE Normalization and Large-Scale Structure			      */
/*      Emory F. Bunn  and  Martin White												      */
/*********************************************************************************/

void COBEnormalize()
{

   const double fourpi=4.0*3.14159265;
	double xlog10,h,xlnh,hc,curv,r;
	double c10,d1,d2,d3,d4,d5,d6,d7,x1,x2,x3,x4,x5,x6,x7;
	double sy,s,sx,sxy,sxx,delt,d1pr,d1ppr;
	
	double d2norm;								// Needed for lensing calculation	
		
   xlog10=log(10.0);
   h=h0/100.0;
   xlnh=log(h);

	// Curvature radius
	if(fabs(omegak)>1.0e-3)
  	{
		hc=2.998e5/h0;                      // h*c. c = 2.338*e5 km/s. h = km/sec/pc 
      curv=-omegak/(hc*hc);					  // Gaussian curvature	
      r=1.0/sqrt(fabs(curv));				  // Radious of curvature	
      }

	// Ensuring scalar to tensor ratio
	if(itflag==1)
		for(int in=1;in<=nn;in++)
			for(int il=2;il<=l[l0];il++)
			{
				printf("rat =  %e\n",rat[il]);
				cltt[il][in]=cltt[il][in]*rat[in];
            clet[il][in]=clet[il][in]*rat[in];
            clbt[il][in]=clbt[il][in]*rat[in];
            clct[il][in]=clct[il][in]*rat[in];
            clkk[il][in]=clkk[il][in]*rat[in];
            cltk[il][in]=cltk[il][in]*rat[in];
            }
	
	// COBE normalization
	// fit the spectrum to a quadratic around C_10 with equal weights in logl

	for(int in=1;in<=nn;in++)
	{
		c10=clts[10][in]+cltt[l[2]][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(clts[l[2]][in]+cltt[l[2]][in])/c10-1.0;
      d2=(clts[l[3]][in]+cltt[l[3]][in])/c10-1.0;
      d3=(clts[l[5]][in]+cltt[l[5]][in])/c10-1.0;
      d4=(clts[l[7]][in]+cltt[l[7]][in])/c10-1.0;
      d5=(clts[l[10]][in]+cltt[l[10]][in])/c10-1.0;
      d6=(clts[l[11]][in]+cltt[l[11]][in])/c10-1.0;
      d7=(clts[l[12]][in]+cltt[l[12]][in])/c10-1.0;
		
      x1=log(1.0*l[2])/xlog10-1.0;
      x2=log(1.0*l[3])/xlog10-1.0;
      x3=log(1.0*l[5])/xlog10-1.0;
      x4=log(1.0*l[7])/xlog10-1.0;
      x5=log(1.0*l[10])/xlog10-1.0;
      x6=log(1.0*l[11])/xlog10-1.0;
      x7=log(1.0*l[12])/xlog10-1.0;

		// Quadratic List Square fit      
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);
      
      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;
		
		printf("%e %e %e\n",delt,d1pr,d1ppr);
		
		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;
		
		// Calculations for lensing
		
		d2norm=c10*(1.1e-9)/clts[1][in]/(2*3.14159265);

		if (ict==2)
		{
			for(int itf=1;itf<=ntf;itf++)
			{
         	d2norm=d2norm*exp(xlnh*4.0);
            anorm[itf][in]=d2norm/((4*3.14159265));
				}
			}	
		 
			
		// C_l normalization; output l(l+1)C_l/twopi
		
		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);		
		for(int il=2;il<=l[l0];il++)
		{
			printf("clts = %e \n",clts[il][in]);
			clts[il][in]=clts[il][in]*c10;
	      cles[il][in]=cles[il][in]*c10;
   	   clcs[il][in]=clcs[il][in]*c10;
   	   cltt[il][in]=cltt[il][in]*c10;
   	   clet[il][in]=clet[il][in]*c10;
   	   clbt[il][in]=clbt[il][in]*c10;
   	   clct[il][in]=clct[il][in]*c10;
   	   clkk[il][in]=clkk[il][in]*c10;
   	   cltk[il][in]=cltk[il][in]*c10;
   	   }
		}
	}
	

	
int main()
{

	double gsnunr;
	double xmnu;																									// Store the neutrino mass
	int irt;

	printf("Enter Choice:\n\t (0) CMB\n\t (1) Transfer Function\n\t (2) Both\n");
	printf("For lensed Cl we need (2)\n\nEnter:   ");
	scanf("%d",&ict);
	
   
	if(ict!=1)																										// If CMB is required then
	{

		printf("\nEnter value of lmax,ketamax (e.g. 1500 3000)");
		printf("\n(Must be consistent with the file)\n(in flat case)");
		printf("\n\nEnter: ");
		scanf("%d %lf",&lmoin,&akmax0);
		
		initlval(lmoin);
		}

	/********************************************************************************************************************************/
	/*	 Output at z=10 saves 50% of time; for CDM models TF remains unchanged between z=0 and z=10; for MDM, lambda, nonflat		  */
	/*	 models one should integrate until z=0 																												  */
	/********************************************************************************************************************************/	

	if(ict!=0)																										// If transfer function is required then	
	{
		printf("\nEnter Transfer Function kmax(h/Mpc),");
		printf(" \n# of k per log. int. (e.g. 5 5)\nEnter: ");
		scanf("%lf %d",&akmaxt,&nlnkt);
		printf("\nEnter the number and the redshifts of the output: (e.g. 1 0)\n");
		printf("If more than one TF is requested then the redshifts");
		printf("\nhave to be specified in the decreasing order\nEnter the redshifts:");
		scanf("%d",&ntf);
		for(int i=1;i<=ntf;i++)
		{
			scanf("%lf",&ztf[i]);
			}
		}
	
	nflag_rho=0;
	printf("\nConstant equation of state with dark energy");
	printf("\n (1) perturbation\n (3) without perturbation\n");
	printf(" (2) table with a,w with perturbation\n (4) table ");
	printf("with a,w without perturbation\nEnter 1 or 3 please : ");
	scanf("%d",&ndyn);		
	
	if(ndyn==1 || ndyn==3)
	{
		printf("\nEnter w_dyn (e.g. -1) :");
		scanf("%lf",&wdyn);
		if(wdyn == -1) ndyn=0;
		}
	else
	{
		printf("Sorry these parts are not coded till now");
		exit(1);
		}
	
	
   printf("\nEnter (1) to include the 5th dimension");
   printf("\nEnter (0) for usual FRW spacetime\nEnter choice :");
   scanf("%d",&dimflag);
   
   printf("\nEnter Omega_b, Omega_c, Omega_de,Omega_nu \n(e.g. 0.05 .3 0.65 0)\n");  		
	scanf("%lf %lf %lf %lf",&omegab,&omegac,&omegav,&omegan);

	omega=omegab+omegac+omegav+omegan;	

	if(dimflag==1) 
	{
		omega=(sqrt(omegav)+sqrt(omegav + omegac + omegab + omegan));
		omega=omega*omega;
		}
	omegak=1.0-omega;

	printf("\nEnter H0, Tcmb, Y_He, N_nu(massless), N_nu(massive), g*(massive)");
	printf("\n (e.g. 72 2.726 0.24 3.04 0 10.75)\nEnter: ");
	scanf("%lf %lf %lf %lf %lf %lf",&h0,&tcmb,&yhe,&annur,&annunr,&gsnunr);

	akmaxt = akmaxt*h0/100;
		
	/********************************************************************************************************************************/			
	/*  Test the consistency of the input values. If not consistent then prompt it   															  */
	/********************************************************************************************************************************/	

   if(h0<25 || h0>100 ) 		 printf("\nWarning: H0 has units of km/s/Mpc. Your value is weird.");
   if(tcmb<2.7 || tcmb>2.8) 	 printf("\nWarning: Tcmb has units of K. Your value is weird.");
   if(yhe<0.2 || yhe>0.3) 		 printf("\nWarning: Y_He is the Helium fraction of baryons. Your value is weird. ");
   if(annur<0 || annur>3.1) 	 printf("\nWarning: N_nu(massless) is strange");
   if(annunr<0 || annunr>3.1)  printf("\nWarning: N_nu(massive) is strange");
   if (annunr<1 && omegan>0.0) printf("\nWarning: N_nu(massive) should be 1, 2, or 3', For non zero omegan\n\n");	
		

	printf("\nEnter (0) for Peebles recombination \n      (1) for recfast\nEnter: ");
	scanf("%d",&rcflag);
	
	if(rcflag==1)
	{
		recfast();
		printf("\nSorry RECFAST recombination function is not present");
		}	

//	FILE *fp33;
//	fp33 = fopen("recdata.d","w");
//	for(int diid=1;diid<=Nz0+1;diid++)
//		fprintf(fp33,"%e\t%e\n",1/(1+zrec[diid]),xrec[diid]);
//	fclose(fp33);
//	exit(1);	

	/**************************************************************************************************************/	
	/*   Evaluate general cosmology constants.																						  */	
	/*																																				  */
	/*   grho gives the contribution to the expansion rate from:																  */
	/*     (g) photons,																														  */	
	/*     (r) one flavor of relativistic neutrino (2 degrees of freedom),													  */
	/*     (m) nonrelativistic matter (for Omega=1).																				  */	
	/*     grho is actually 8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).													  */	
	/*     a=tau(Mpc)*adotrad, with a=1 today																							  */			
	/*     (Used only to set the initial conformal time.)																			  */	
	/**************************************************************************************************************/	
	
	grhom=(3.3379e-11)*h0*h0;                                           // 3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
	grhog=(1.4952e-13)*pow(tcmb,4);                                     // 8*pi*G/c^2*4*sigma_B/c^3 T^4        
	grhor=(3.3957e-14)*pow(tcmb,4);                                     // 7/8*(4/11)^(4/3)*grhog (per neutrino species)


	if(annunr!=0)
	{
		grhonr=grhog*0.875*pow((43/11)/gsnunr,(4/3));
		xmnu=omegan*(h0*h0/10000.0)*(gsnunr/10.75)*(93.5/annunr);
		printf("Nutrino mass = %e ev",xmnu);													// Calculating Neutrino mass
		}
	
	adotrad = sqrt((grhog+grhor*annur)/3);													// adotrad gives the relation a(tau) in the radiation era
	
															

	printf("\n\nEnter  (0) for no reionization \n");
	printf("       (1) for specified optical depth to lss(xe=1)\n");
	printf("       (2) for specified redshift and xe\n");
	printf("\nEnter Choice: ");
	scanf("%d",&riflag);
	
   zri=0.0;
   rif=0.0;
   optdlss=0.0;	

   if(riflag==1)
   {
      printf("Enter optical depth to lss\n");
      scanf("%lf",&optdlss);
      rif=1.0;
      }

   if (riflag==2)
   {
      printf("\nEnter redshift, ionization fraction (e.g. 50 0.2) : ");
      scanf("%lf %lf",&zri,&rif);
      zristp=0.07*zri-1.0;
      if(zristp<0.0) zristp=0.0;
      }
   
      if (ict!=1)
   { 
      printf("\nEnter :\n (0) for scalar alone \n");
      printf(" (1) for tensor+scalar\n (2)");
      printf(" for tensors alone\n (3) for scalar ");
      printf("(k<k*)\n (4) for scalar (k>k*)\nEnter choice: ");
      scanf("%d",&itflag);
    	}

	/********************************************************************************************************************************/   
	/* K split																																							  */	
	/********************************************************************************************************************************/
	aksplit=1.0;
	kcutflag=0;
	
	if(itflag==3)
	{
		itflag=0;
		kcutflag=1;
		}
	
	if(itflag==4)
	{
		itflag=0;
		kcutflag=-1;
		}
	
	printf("\nEnter the inflationary model :\n");
	printf("  (0)  Single field inflation\n  (1)  Two field inflation\n");
	scanf("%d",&twofield);	

	if(twofield == 0)
	{
		if(itflag!=2)
		{
			printf("Number and values of scaler spectral index n,");
			printf("\nand its running alpha_n (e.g 1 1 0)\nEnter choice: ");
			scanf("%d",&nn);
			for(int i=1;i<=nn;i++)
			{
				scanf("%lf %lf",&an[i],&alphans[i]);
				}
			}	
		
		if(itflag==1)
		{
			int itn;
			printf("\nTensor spectral index given by \n (0)");
			printf(" nt=ns-1  \n (1) different\nEnter choice: ");
			scanf("%d",&itn);
			if(itn==0)
				for(int i=1;i<=nn;i++)
				{
					ant[i]=an[i]-1;
					alphant[i]=alphans[i];
					}
			else
			{
				printf("Values of tensor spectral index n,");
				printf("\n and its running alpha_n\n (e.g 0 0)");
				for(int i=1;i<=nn;i++)
				{
					scanf("%lf %lf",&ant[i],&alphant[i]);
					}	
				}
			}

		if(itflag==2)
		{
			printf("\n\nNumber and values of scaler spectral");
			printf(" index n,\n and its running alpha_n\n (e.g 1 0 0)");
			scanf("%d",&nn);
			for(int i=1;i<=nn;i++)
			{
				scanf("%lf %lf",&ant[i],&alphant[i]);
				}
			}
			
		if(itflag==1)
		{
			printf("\nRatio of tensor to scalar quadrupole");
			printf(" given by\n  (0) 7*(1-n_S)  \n  (1) different\nEnter: ");
			scanf("%d",&irt);
			if(irt==0)
				for(int i=1;i<=nn;i++)
				{
					rat[i]=7.0*(1-an[i]);
					}
			else
			{
				printf("\nEnter the values of T/S:");
				for(int i=1;i<=nn;i++)
					scanf("%lf",&rat[i]);
				}		
		
			}
	
		if(itflag==0)
		{
			for(int i=1;i<=nn;i++)
				rat[i]=1.0;
			}
	
	
		lensflag=0; 
		
		if(itflag!=2 && ict==2)
		{
			printf("\nEnter\n (0) for unlensed Cl only");
			printf("\n (1) for lensed Cl linear evalution\n");
			printf(" (2) for lensed Cl nonlinear evalution");
			if(itflag==1)
				printf("\nOnly Scaler Cl's are lensed");
			printf("\nEnter Choice: ");
			scanf("%d",&lensflag);
			}		
			
		if(lensflag!=0 && ict==0)
		{
			printf("You did not request the transfer function");
			printf(" and so calculation of lensing can not be");
			printf(" done. Restart:");
			exit(1);
			}	
	
		if(itflag!=2)
		{
			printf("\nEnter initial conditions\n"); 
			printf("	(1) Isentropic (adiabatic)\n");
			printf(" (2) Isocurvature CDM\n");
			printf(" (3) Isocurvature baryon\n");
			printf(" (4) Isocurvature seed conditions");
			printf("\nEnter choice: ");
			scanf("%d",&initfl);
			}
		else
			initfl=1;
		}
	
	if(twofield == 1)
	{
		initfl = 1;
		nn = 1;
	
		printf("\nDo you want the standard model inplimented in the code. ");
		printf("\nOr you want to input the power from some file :\n");
		printf(" (0) Use function of the code \n");
		printf(" (1) Use data from file\n");
		scanf("%d",&powerfileflag);
		if(powerfileflag == 1)
			readpowers();
		if(powerfileflag == 0)
		{
			printf("Enter sH,s0,R \n       (60.0 50.0 5.0)\n");
			scanf("%lf %lf %lf",&powerssH,&powerss0,&powersR);
			powersh00 = h0;
			}
		}

	printf("akmax0 = %e",akmax0);
	initjl();		printf("done 1");
//	exit(1);

	cmbflat();
	printf("done");
	if(twofield == 0)
		COBEnormalize();
	else if(twofield == 1)	
		COBEnormalizetwofield();

	if(lensflag!=0)
		lensing();
	
	FILE *sanfp,*sanfpt,*sanfpe,*sanfpc;
	sanfp = fopen("cl.d","w");	
	sanfpt = fopen("clt.d","w");	
	sanfpe = fopen("cle.d","w");	
	sanfpc = fopen("clc.d","w");	

	if(twofield == 0)
		for(int in=1;in<=nn;in++)
			for(int il=2;il<=l[l0];il++)
				fprintf(sanfp,"%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n",clts[il][in],cles[il][in],clcs[il][in],clkk[il][in],cltk[il][in],cltt[il][in],clet[il][in],clbt[il][in],clct[il][in]);	

	else if(twofield == 1)
		for(int in=1;in<=nn;in++)
			for(int il=2;il<=l[l0];il++)
			{
				fprintf(sanfpt,"%e\t %e\t %e\t %e\t\n",cltsadb[il][in],cltsiso[il][in],cltscross[il][in],cltstotal[il][in]);
				fprintf(sanfpe,"%e\t %e\t %e\t %e\t\n",clesadb[il][in],clesiso[il][in],clescross[il][in],clestotal[il][in]);
				fprintf(sanfpc,"%e\t %e\t %e\t %e\t\n",clcsadb[il][in],clcsiso[il][in],clcscross[il][in],clcstotal[il][in]);	
				}
	else 
		printf("\nWrong choice to inflation model");


	fclose(sanfp);
	printf("Finally done");
	}