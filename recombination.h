#ifndef	_RECOMBINATION_H_
#define	_RECOMBINATION_H_

#include "variables.h"
#include "numericx.h"
#include "initval.h"
#include "neutrino.h"
#include "darkenergy.h"

double tauminn,dlntau;														// tauminn <- minimum conformal time from which we will start the calculation
																					// dlntau <- time step in the log scale	... For calculating the recombination (before greading)

/*************************************************************************************************************************************/
/*  Integrate the ionization fraction xe for hydrogen semi-implicitly from tau to tau+dtau, treating tempb,								  	 */
/*   a, and adot as constants during this interval.																											 */
/*************************************************************************************************************************************/

double ionize(double tempb,double a,double adot,double dtau,double xe)
{
																					
																					//  Switch for implicit (switch=1.0) or semi-implicit (switch=0.5) scheme.
   double switch1=0.5;														//  Recombination coefficient (in sqrt(K)/Mpc).
   double alpha0=(2.3866e-6)*(1-yhe)*omegab*h0*h0;					//  Coefficient for correction of radiative decay (dimensionless).
 	double crec,phi2;
  	double alpha,beta,cpeebles;
  	double cp1,cp2,bb,b1,bbxe,aa;
  	double rat;


   double dec2g=8.4678e14;													//  Two-photon decay rate (in 1/Mpc)
   double tion=1.5789e5,beta0=43.082;									//  Ionization temperature and coefficient.
   
   phi2=0.448*log(tion/tempb);
   phi2=max(phi2,0.0);
   
	/*************************************************************************************************************/
	// Variables Details
	//
	// dec2g  =  Lambda_2s,1s 
	//        =  8.227 sec^-1 =  8.227 * 1.029272e14 Mpc^-1 = 8.4678e14 Mpc^-1
	//
	// tion = B1/kb = 13.6 ev * 1.6*10^-19 / 1.38 * 10^-23 = 1.5789e5 Kelvin
	//
	//	phi2 = 0.448*log(B1/(kb*Tb))										// eq(79) Ma & Bertschinger
	//
	//
	// Lalpha =  1.215670018e-07 m;
	// K = Lalpha^3/(8.0*pi);
	// Mpc = 3.085678e22 m;
	// m_H = 1.673575e-27 kg;
	// G=6.67428e-11 m^3/kg/s^2;
	// K*3*(1.0e3/Mpc)^2/(8*pi*G*m_H) = 8.0230194e-26				// The square term came from the Hubble constant
	// 
	// Constant = (m_e*k_b/2/pi/h_bar)^1.5 * ((64*pi/(27*pi)^0.5)*(e^4/(m_e^2 * c^3))) * (k_b/B1)^-0.5 * 1/(Mpc2sec)
	// B1 = m_e * e^4 / h_bar^2
	// => Constant = 5.13*10^18 = exp(43.082)
	// beta0 <- 43.082														// beta0 is logarithm value of the multiplication constant 
	//																				// eq(78) & eq(79) Ma & Bertschinger
	//
	//
	// HO = 1.0D3/Mpc;
	// c = 2.99792458e8 cm/sec;
	// MPC_in_sec = Mpc/c;
	// e = 4.803e-10 esu;
	// me = 9.10938215e-28 gm;
	// G = 6.674e-8 cm^3 /gm/s^2;
	// mh = 1.673575e-24 gm;
	// (3/(8*pi*G*mh))*((64*pi/sqrt(27*pi))*(e^4/(me^2*c^3))*sqrt(13.6/(8.617343e-5)))*(HO)^2*MPC_in_sec = 2.3848e-006 /sec
	//
	/***************************************************************************************************************/
  		
   crec=(8.023019e-26)*(1-yhe)*omegab*h0*h0;							// +eq(81) Ma & Bertschinger 
   																				// 8.0138e-26 was given in CMBFAST

   alpha=(alpha0/sqrt(tempb))*(phi2/(a*a*a));
   beta=tempb*phi2*exp(beta0-tion/tempb);								//  Peebles' correction factor.
   
   if (tempb<200.0)
   {
   	cpeebles=1.0;
   	}
   else
   {
   	cp1=crec*dec2g*(1.0-xe)/(a*adot);
      cp2=crec*tempb*phi2*exp(beta0-0.25*tion/tempb)*(1.0-xe)/(a*adot);
      cpeebles=(1.0+cp1)/(1.0+cp1+cp2);
   	}
																					//  Integrate dxe=bb*(1-xe)-aa*xe*xe by averaging rhs at current tau
																					//  (fraction 1-switch) and future tau (fraction switch).
	aa=a*dtau*alpha*cpeebles;
   bb=a*dtau*beta*cpeebles;
   
   /***************************************************************************************************************/
   //  Integration method
   //  f(t1) = a*(1-x1) + b*x1^2
   //  (x2-x1)/dt = (1-s)*f(t2)+s*f(t1)
   //  Solve for x2 using Sridharacharya rule in terms of x1
   /***************************************************************************************************************/
   
   
   b1=1.0+switch1*bb;
   bbxe=bb+xe-(1.0-switch1)*(bb*xe+aa*xe*xe);
   rat=switch1*aa*bbxe/(b1*b1);											//  Prevent roundoff error.
        
   if (rat < 1.0e-6)
   	xe=(bbxe/b1)*(1.0-rat);
   else
   	xe=b1/(2.0*switch1*aa)*(sqrt(4.0*rat+1.0)-1.0);

	return(xe);
	}        
	
/************************************************************************************************************************************/
//	Saha recombination equation
// Here x_H0 will return ionization fraction of the Hydrogen atom
/************************************************************************************************************************************/
double ionsaha(double tempb,double a,double xe) 
{
   double mu_H = 1.0/(1.0-yhe);     	//Mass per H atom
	double HO = h0/(3.0856775807e19);  	// in seconds (that's originals recfast's HO)
  
  //  today's number density of hydrogen
   double Nnow = 3*HO*HO*omegab/(8.0*3.141592653589*G*mu_H*m_H);
	double CB1 = h_P*C*L_H_ion/k_B;
   double CR = 2.0*3.141592653589*(m_e/h_P)*(k_B/h_P);
   double rhs = exp( 1.5 * log(CR*tempb*a*a) - CB1/(tempb) ) / Nnow;
   double x_H0 = 0.5 * (sqrt(rhs*rhs +4.0*rhs ) - rhs );
	if(x_H0<0.0001) x_H0 +=1e-5;
   return x_H0;
   }
   
   
/************************************************************************************************************************************/
/*  Compute the helium ionization fractions using the Saha equation.																						*/
/*  x0 is the hydrogen ionization fraction n(H+)/n(H) (input),																								*/
/*  x1 is the helium first ionization fraction n(He+)/n(He)	(input and output), and																	*/
/*  x2 is the helium second ionization fraction n(He++)/n(He) (input and output).																	*/
/************************************************************************************************************************************/

void ionhe(double tempb,double a,double xsan[])
{
	double tion1=2.855e5,tion2=6.313e5;								// 1st and 2nd Ionization temperature of He
	double b0=(2.150e24)/((1.0-yhe)*omegab*h0*h0);  			//  Constant for electron partition function per baryon.
	
	double xe,x1,x2,x0,r1,r2,b;
	
	double x2new,error=1.0;
	double cons;
	
	x0=xsan[0];
	x1=xsan[1];
	x2=xsan[2];
	
	/********************************************************************************************************************************/	
	// Mpc = 3.085678e22;
	// m_H = 1.673575e-24 gm;
	// G=6.67428e-8 cm^3 /gm/s^2;
	// denom = (3*(1.0e3/Mpc)^2)/(8*pi*G*m_H)
	// me = 9.10938215e-28 gm;
	// kb = 1.380650e-16 erg/K;
	// hb = 1.054571628e-27 erg/sec;
	// 
	// ((me*kb/(2*pi*hb*hb))^1.5)/denom = 2.1514e24				// eq(76) Ma & Bertschinger
	/********************************************************************************************************************************/	
	
																				// eq(76) Ma & Bertschinger
   b=b0*a*a*a*tempb*sqrt(tempb);										//  Electron partition function per baryon.
																				//  Dimensionless right-hand sides in Saha equations.
   if (fabs(tion1/tempb)<100.0) 
      r1=4.0*b*exp(-tion1/tempb);
   else
      r1=0.0;

   if (fabs(tion2/tempb)<150.0)
   {
	   r2=b*exp(-tion2/tempb);
      }			
   else
   {
	   r2=0.0;        
      }																		//  Solve coupled equations iteratively.
   cons=0.25*yhe/(1.0-yhe);

	/********************************************************************************************************************************/
	// Numerical method involved in solving the equation
	//
	// e1:   x2/x1 = r2/xe
	// e2:   x1/(1-x1-x2) = r1/xe                            // eq(76) Ma & Bertschinger
	//
	// (e2) will give 
	// x1 = xe*r1/(r1*r2+xe*r1+xe*xe)
	// So x2 will be 
	// x2new=r1*r2/(r1*r2+xe*r1+xe*xe)
	// Itrerate this till it converge to some fixed value									
	/********************************************************************************************************************************/
	
	while(fabs(error)>(1.0e-12))
	{
  		xe=x0+cons*(x1+2*x2);
   	x2new=r1*r2/(r1*r2+xe*r1+xe*xe);
		x1=xe*r1/(r1*r2+xe*r1+xe*xe);
		error=x2-x2new;
		x2=x2new;
	  	}
   xsan[0]=x0;
   xsan[1]=x1;
   xsan[2]=x2;	
   }

/*********************************************************************************************************************************/
/*     Compute and save unperturbed baryon temperature and ionization fraction as a function of time.  With nthermo=10000,			*/
/*     xe(tau) has a relative accuracy (numerical integration precision) better than 1.e-5.													*/
/*********************************************************************************************************************************/

int finithermo(double taumin,double taumax,double dlntau0,int nstepsan[],double tausann[])
{
	double tau01;
	double rxsan[2];
	double tb[nthermo+1];
	double tbhalf;
	double thomc,fe,etc,a2t;
	double rhonu;																					// Neutrino density and pressure
	double grho,d;
	int j2ri3,j2ri2;																				// 
	double dtri,dtri0,dtauri,dtauri0;
	double atemp[nthermo+1],datemp[nthermo+1],xe[nthermo+1];											
	double a2,adot,tg0;
	double vfi,cf1;
	double ahalf,adothalf;
	double tmp1,tmp2;
	int iv,ns,ncount;
	double tgh,xod;
	double accstep,fact;
	double a,a0,adot0;
	double tautemp,taurist1,taurend2,taurend1;
	double vismax,tvis;																				// Visibility function
	double dtau,tau,xe0,x1,x2;
	double xsan[3];
	double a02;
	double ddopac;
	int nstep;
	double taurend;		
	double dtbdla;																						// We have to return taurend
	
	
	int i;																								// Will be used by for loops and from some other places.. tempurary variable
	int n1;																								// Number of gridpoints during recombination	
	
	double akthom;																						// Np_now * sigma_thomson ... for calculating optical depth 	
	double thomc0;																						// Used for calculating the baryon temperature
	double barssc0=9.1843e-14, barssc;															// Will be used for calculating baryon sound speed
																											// barssc0 = k_B / m_p / c^2 = 9.1843e-014
																											//           CMBFAST uses 9.1820e-014
	FILE *fp1;
	
	fp1=fopen("resaha.d","w");
	
	nstep = nstepsan[0];
	ncount=0;
	
	
	// fp1: Initiallizing variables
	
   thomc0=(5.0577e-8)*pow(tcmb,4);
   akthom=(2.3048e-9)*(1-yhe)*omegab*h0*h0;													// Np_now * sigma_thomson

	/********************************************************************************************************************************/   
   //
   // akthom = Np_now * sigma_thomson
   // 
   // Mpc = 3.085678e22;
	// m_H = 1.673575e-24;
	// G=6.67428e-8;
	// sigma_thomson = 6.6524616e-25;
	// (3*(1.0e3/Mpc)^2)/(8*pi*G*m_H)*100*Mpc*sigma_thomson = 2.3048e-9
   //
   //
   // thomc0 = MPC_in_sec*(8/3)*(sigma_thomson/(m_e*c))*a  * T^4                  // a is the radiation constant for u = a.T^4 (Not to be confused with scale factor)
   //        = 5.0577e-8 * Tcmb^4
   //			 ..... it will be used for calculating the baryon temperature
	/********************************************************************************************************************************/   


	// fp2: Initial conditions: assume radiation-dominated universe.
	tauminn=0.05*taumin;
	adot0=adotrad;
   a0=adotrad*tauminn;
   a02=a0*a0;
   
   tb[1]=tcmb/a0;	
																											// Assume that any entropy generation occurs before tauminn.
   xe0=1.0;																								// This gives wrong temperature before pair annihilation, but the error is harmless.
   x1=0.0;
   x2=1.0;
   xe[1]=xe0+0.25*yhe/(1.0-yhe)*(x1+2*x2);
   barssc=barssc0*(1.0-0.75*yhe+(1.0-yhe)*xe[1]);
   cs2[1]=(4.0/3.0)*barssc*tb[1];
   dotmu[1]=xe[1]*akthom/a02;
   sdotmu[1]=0;
   atemp[1]=a0;
   dlxedla[1]=0.0;
		
	tau01=tauminn;
	dlntau=log(tau0/tauminn)/(nthermo-1);
	
	//fp3: Calculating the ionization fraction
	
	for(int i=2;i<=nthermo;i++)
	{
		tau=tauminn*exp((i-1)*dlntau);
		dtau=tau-tau01;

		/***********************************************************************************************************/
		//             Integrate Friedmann equation using inverse trapezoidal rule. (PC method)						  //
		/***********************************************************************************************************/
		
      a=a0+adot0*dtau;																		
      
      a2=a*a;
      
      nu1(a,rxsan);
		rhonu=rxsan[0];																		// Calculate neutrino density	
	
      omegavdyn = omegav*dynrho(a);

      if(dimflag==0)
      	grho=grhom*(omegac+omegab)/a+(grhog+grhor*annur+grhonr*annunr*rhonu)/a2 + grhom*omegavdyn*a2;
		else	
		{
      	grho=(sqrt(grhom*(omegac+omegab)/a+(grhog+grhor*(annur+annunr*rhonu))/a2+grhom*omegav*a2)+sqrt(omegav*grhom)*a);
      	grho=grho*grho;
			}
					
		adot=sqrt(grho/3.0)*a;
      a=a0+2.0*dtau/(1.0/adot0+1.0/adot);
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
      atemp[i]=a;
		
		tg0=tcmb/a0;
      ahalf=0.5*(a0+a);
      adothalf=0.5*(adot0+adot);
		 
      fe=(1.0-yhe)*xe[i-1]/(1.0-0.75*yhe+(1.0-yhe)*xe[i-1]);							// fe=number of free electrons divided by total number of free baryon + electron
      thomc=thomc0*fe/(adothalf*pow(ahalf,3));												// particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
      etc=exp(-thomc*(a-a0));																		// more accuracy is required (unlikely) then this can be iterated with
      a2t=a0*a0*(tb[i-1]-tg0)*etc-tcmb/thomc*(1.0-etc);									// the solution of the ionization equation.
      tb[i]=tcmb/a+a2t/(a*a);
      																									// Integrate ionization equation.
      tbhalf=0.50*(tb[i-1]+tb[i]);

      if ((zri!=0.0)&&(tau>(9.0*taurist/10.0))) 											// If there is re-ionization, smoothly increase xe to the requested value.
      {  

      	if(ncount==0)
         	ncount=i-1;
         	
	      xod=150.0*(tau-taurist)/taurist;														//	SMOOTH REIONIZATION
				        	
         if (xod>90)
         	tgh=1.0;
         else
         	tgh=(exp(xod)-exp(-xod))/(exp(xod)+exp(-xod));
         
         xe[i]=(rif-xe[ncount])*(tgh+1.0)/2.0+xe[ncount];
         dlxedla[i]=(xe[i]-xe[i-1])*a/(dtau*adot*xe[i]);
         }
         
      else
      {
      	if(rcflag==0 || rcflag==2)
      	{
				if(rcflag==0)
	          	xe0=ionize(tbhalf,ahalf,adothalf,dtau,xe0);
   			else
	          	xe0=ionsaha(tbhalf,ahalf,xe0);
	      	xsan[0]=xe0;
	      	xsan[1]=x1;
	      	xsan[2]=x2;
          	ionhe(tb[i],a,xsan);

            xe0=xsan[0];
            x1=xsan[1];
            x2=xsan[2];
            xe[i]=xe0+0.25*yhe*(x1+2*x2)/(1.0-yhe);
            dlxedla[i]=(((xe[i]-xe[i-1])/xe[i])/dtau)*(a/adot);
            }
         else
         {
	         printf("This part is not done till now. Try Peebles recombination.");
         	//call recint(a,xe[i]);
            //dlxedla[i]=(xe[i]-xe[i-1])*a/(dtau*adot*xe[i]);
				}
         }
	   
	        
      barssc=barssc0*(1.0-0.75*yhe+(1.0-yhe)*xe[i]);
      dtbdla=-2.0*tb[i]-thomc*adothalf/adot*(a*tb[i]-tcmb);
      cs2[i]=(1.0 + dtbdla/tb[i]/3.0)*barssc*tb[i];															// Sound speed square
      dotmu[i]=xe[i]*akthom/a2;																						// Calculation of the d_mu/d_tau
		
	   fprintf(fp1,"%e\t%.8e\t%.8e\t%.8e\t%.8e\n",a,xe[i],xe0,x1,x2);
      a0=a;
      tau01=tau;
      adot0=adot;
		}
	fclose(fp1);
//	printf("San");
//	exit(1);
	
	if((xe[nthermo]<rif)&&(zri!=0.0))
	{
		printf("\nWarning: We use a smooth function to approach your specified reionization");
      printf("\nfraction. The redshift that is deduced from your input parameters is so low that our");
      printf("\nsmooth function does not reach the required value. You should go in to subroutine finithermo");
      printf("\nand play with the shape of this smooth function. Search for SMOOTH REIONIZATION for");
      printf("\nthe place where the function is set.");
      }

																											// Saving tau values and the quantities needed to calculate
																											// the derivatives of the visibility function appearing in the sources.
   for(int j1=nthermo-1;j1>=1;j1--)
   {
	   tmp1=dotmu[j1]*tauminn*exp((j1-1)*dlntau);
      tmp2=dotmu[j1+1]*tauminn*exp((j1)*dlntau);											// Here the integration is done by (d_mu/d_tau)*tau*(d_tau/tau)... p.s. dlogtau = d_ tau/tau
      sdotmu[j1]=sdotmu[j1+1]+0.5*(tmp1+tmp2)*dlntau;										// Sdotmu = mu <- Optical depth at conformal time tau .. S stands for integration
      emmu[j1]=exp(-sdotmu[j1]);																	// emmu <- exp(-mu)
      }
	
	emmu[nthermo] = 0;

	// fp4: Fixing the points for greading  
	
   iv=0;
   vfi=0.0;
																											// Getting the starting and finishing times for decoupling.
   if (ncount==0) 
   {
	   cf1=1.0;
      ns=nthermo;
      }
   
   else
   {
	   cf1=exp(sdotmu[nthermo]-sdotmu[ncount]);
      ns=ncount;
      }
   

   for(int j1=1;j1<=ns;j1++)
   {
	   tau=tauminn*exp((j1-1)*dlntau);
      vfi=vfi+emmu[j1]*dotmu[j1]*cf1*dlntau*tau;
      if ((iv==0)&&(vfi > 1.0e-12))
      {
	      taurst=(9.0/10.0)*tauminn*exp((j1-1)*dlntau);
       
         iv=1;
         }
         
      if ((iv==1)&&(vfi > 0.999))
      {
	      taurend1=(tauminn*exp((j1-1)*dlntau));
         taurend1=max(taurend1,taurend1*sqrt(2500.0/((omegac+omegab)*(h0*h0)))); // Most probably this step is not needed. It will only slow down the code
         iv=2;
         }
		}
   if (iv!=2)	
   {
	   printf("\n Failed to find the end of the recombination process.");
	   printf("\n Recombinaion end is taken as 1.3 times start of reionization");	// Calculating the timesteps during recombination.
   	taurend1=1.5*(tauminn*exp((ncount-1)*dlntau));										// A small timestep is needed because of the low ks, the cancellation
		}																									// across the visibility function are not good otherwise
 																											// In CMBFAST this number is taken as 2.5*(tauminn*exp((ncount-1)*dlntau));
   accstep=3.0;

	if (dtaurec!=0.0)
	{
		dtaurec=min(dtaurec,taurst/(40.0*accstep));
		}
	else
	{
		dtaurec=taurst/(40.0*accstep);
      }
      
      
	taurend2=dtaurec/dlntau0;																		// For smooth greading... At the end of recombination dlntau0 = dtaurec/taurec
									
																											// In models where reionization starts very early so that it cut into what we call
   taurend=max(taurend1,taurend2);																// recombination we will make the timesteps there at least as small as the ones 
   taurend=min(taurend,9.0*taurist/10.0);														// we were using during recombination. If not timesteps 1.5 times bigger suffice.

   fact=1.5*accstep;
   if (taurend==(9.0*taurist/10.0))
   	fact=1.0*accstep;
 	
 	// Modify (dtaurec)	..... timestep during recombination
   n1=(int)((taurend-taurst)/dtaurec+1);
   dtaurec=(taurend-taurst)/(n1-1);																// Calculating the timesteps after recombination (logarithmic outside  
   n1=n1+1;																								// re-ionization scattering surface).
	
 	// Modify (dlntau0)	..... timestep between recombination and reionization
	nstep=(int)(log(taumax/taurend)/dlntau0);
   dlntau0=log(taumax/taurend)/nstep;
   nstep=nstep+n1;
  

   nri0=50;																						 		// Adjusting if there is reionization. There will be nri0 points to sample th//e
																											// quick rise in the free electron density. After that, timesteps of length  
   if (zri!=0)																							// dtauri until tauristp.
   {
		taurist1=taurist;
      j2ri1=(int)(log(9.0*taurist1/(taurend*10.0))/dlntau0+n1);
	
		tautemp=min(21.0*taurist1/20.0,tauristp);
      j2ri2=(int)(log(tautemp/taurend)/dlntau0+n1);
	
		j2ri3=(int)(log(tauristp/taurend)/dlntau0+n1);
      dtri0=taurend*(exp(dlntau0*(j2ri2-n1))-exp(dlntau0*(j2ri1-n1)));
      dtri=taurend*(exp(dlntau0*(j2ri3-n1))-exp(dlntau0*(j2ri2-n1)));
      dtauri0=dtri0/nri0;
	
	
      dtauri=dtaurec*fact;
      dtauri=min(dtauri,dtri/10.0);
	
		if (dtauri > 0.0)
		{
			nri=(int)(dtri/dtauri)+1;
         dtauri=dtri/nri;
         }
      else
      {
	      nri=0;
         dtauri=0.0;
         }
      nstep=nstep+nri0+nri+j2ri1-j2ri3;
      }
	
	else
	{
		j2ri1=0;
      j2ri2=0;
      }
	if (nstep > nstep0)
	{
		printf("Sorry, the arrays were dimensioned for a maximum of %d, timesteps.",nstep0);
		printf("The model you requested needs %d, Please make the arrays bigger by ",nstep);
		printf("making nstep0 bigger where it appears");
		}
	
	// fp5: Calculating derivative of the values stored in the arrys for interpolating them in the new grid points 
	
	splini();
   splder(cs2,dcs2,nthermo);
   splder(dotmu,ddotmu,nthermo);
   splder(ddotmu,dddotmu,nthermo);
   splder(dddotmu,ddddotmu,nthermo);
   splder(emmu,demmu,nthermo);
   splder(dlxedla,ddlxedla,nthermo);
   splder(atemp,datemp,nthermo);

   vismax=0.0;																						// Find the peak of the visibility function directly in this arrays:
   
   for(i=1;i<=nthermo;i++)
   {
	   tvis=dotmu[i]*emmu[i];
	   
	   if (tvis > vismax)
	   {
		   vismax=tvis;
         tau=tauminn*exp((i-1)*dlntau);
         taurmax=tau;
         armax=atemp[i];
         }
		}
	
   if(armax > 0.003)
   {
	   printf("The maximum of the visibility function occurs at z=%e",(1/armax-1));
      printf(" which looks odd. If you are using the k splitting technique and");
      printf("you are computing a model with very high tau you should worry about the accuracy");
      }


   vismax=0.0;

   for(int j2=2;j2<=nstep;j2++)
   {
		// fp6: Calculate the conformal time 

	   if (j2<=n1)
	   {
         tau=taurst+(j2-2)*dtaurec;
      	}
      else
      {
         if ((zri==0)||(j2<=j2ri1))
         {
            tau=taurend*exp(dlntau0*(j2-n1));
         	}
         else
         {
            if (j2 < (j2ri1+nri+nri0))
            {
               if (j2 <= (j2ri1+nri0))
               {
                  tau=atau0[j2ri1]+dtauri0*(j2-j2ri1);
               	}
               else
               {
                  tau=atau0[j2ri1+nri0]+dtauri*(j2-j2ri1-nri0);
             		}
            	}	
            else
            {
               tau=taurend*exp(dlntau0*(j2-j2ri1-nri0-nri+j2ri3-n1));
           		}
				}   
      	}

       	
      atau0[j2]=tau;
		   
      d=log(tau/tauminn)/dlntau+1.0;															// Cubic-spline interpolation.
      i=(int)d;
      d=d-i;
      
		// fp7: Calculate interpolated values in the grid points  
      
      if (i<nthermo)
      {
      	opac[j2]=dotmu[i]+d*(ddotmu[i]+d*(3.0*(dotmu[i+1]-dotmu[i])-2.0*ddotmu[i]-ddotmu[i+1]+d*(ddotmu[i]+ddotmu[i+1]+2.0*(dotmu[i]-dotmu[i+1]))));
         dopac[j2]=(ddotmu[i]+d*(dddotmu[i]+d*(3.0*(ddotmu[i+1]-ddotmu[i])-2.0*dddotmu[i]-dddotmu[i+1]+d*(dddotmu[i]+dddotmu[i+1]+2.0*(ddotmu[i] -ddotmu[i+1])))))/(tau*dlntau);
         ddopac=(dddotmu[i]+d*(ddddotmu[i]+d*(3.0*(dddotmu[i+1]-dddotmu[i])-2.0*ddddotmu[i]-ddddotmu[i+1]+d*(ddddotmu[i]+ddddotmu[i+1] +2.0*(dddotmu[i]-dddotmu[i+1]))))-(dlntau*dlntau)*tau*dopac[j2])/pow((tau*dlntau),2);
         expmmu[j2]=emmu[i]+d*(demmu[i]+d*(3.0*(emmu[i+1]-emmu[i])-2.0*demmu[i]-demmu[i+1]+d*(demmu[i]+demmu[i+1]+2.0*(emmu[i]-emmu[i+1]))));

         vis[j2]=opac[j2]*expmmu[j2];
																											// Find the maximum of the visibility function in the grid of
																											// used to compute the sources. This is not used anywhere other
																											// than in testing outputs. The values used for shifting the
																											// spectra were computed using the full vectors of lenght nthermo.	
     	 	if (vis[j2]<vismax)
      	{
         	vismax=vis[j2];
        		jrmax=j2;
         	}

      	dvis[j2]=expmmu[j2]*(pow(opac[j2],2)+dopac[j2]);
   	   ddvis[j2]=expmmu[j2]*(pow(opac[j2],3)+3*opac[j2]*dopac[j2]+ddopac);
   		}
   	else
	   {
      	opac[j2]=dotmu[nthermo];
   	   dopac[j2]=ddotmu[nthermo];
 	      ddopac=dddotmu[nthermo];
      	expmmu[j2]=emmu[nthermo];
   	   vis[j2]=opac[j2]*expmmu[j2];
	      dvis[j2]=expmmu[j2]*(pow(opac[j2],2)+dopac[j2]);
      	ddvis[j2]=expmmu[j2]*(pow(opac[j2],3)+3.0*opac[j2]*dopac[j2]+ddopac);
			}

		}
		 
   atau0[1]=0.0;

  
   atau0[nstep]=min(tau0,atau0[nstep]);
																									// saving the length of the timesteps for the time integration.
	// fp7: Calculate the dtau at each point 
   for(int j2=2;j2<=nstep-1;j2++)
   {
      dtau1[j2]=atau0[j2+1]-atau0[j2];
      dtau2[j2]=fabs(atau0[j2+1]-atau0[j2-1])/2.0;
 		}
   dtau1[nstep]=dtau1[nstep-1]/2.0;
   dtau2[nstep]=(atau0[nstep]-atau0[nstep-1])/2.0;
   dtau1[1]=(atau0[2]-atau0[1])/2.0;
   dtau2[1]=dtau2[2]/2.0;

  	tausann[0]=taurend;
  	
   nstepsan[0] = nstep;	
  		
   return(n1);
	}


/************************************************************************************************************************************/
/* 	Compute unperturbed baryon temperature, sound speed squared,																						*/
/*	 	and ionization fraction by interpolating pre-computed tables.																						*/
/************************************************************************************************************************************/
//	cs2b=san[0];
//	opac=san[1];
//	dxe=san[2];
	
void thermo(double tau,double san[])
{
	double cs2b,opac2,dxe,d2;
	int i;
	cs2b=san[0];
	opac2=san[1];
	dxe=san[2];
	
	d2=log(tau/tauminn)/dlntau+1.0;
   i=int(d2);																							
   d2=d2-i;
   
//   for(int j=1;j<=nthermo;j++)
//   {
//	   printf("%d cs %e opac %e dlxedla %e \n",j,dcs2[j],ddotmu[j],ddlxedla[j]);
//	   }
//   exit(1);
   
   if (i<1) 																	// Linear interpolation if out of bounds (should not occur).
   {	
	   cs2b=cs2[1]+(d2-1)*dcs2[1];
      opac2=dotmu[1]+(d2-1)*ddotmu[1];
      dxe=dlxedla[1]+(d2-1)*ddlxedla[1];
      }
		
	else if (i>=nthermo)
 	{	
	 	cs2b=cs2[nthermo]+(d2-nthermo)*dcs2[nthermo];
      opac2=dotmu[nthermo]+(d2-nthermo)*ddotmu[nthermo];
      dxe=dlxedla[nthermo]+(d2-nthermo)*ddlxedla[nthermo];
      } 	
	else																			//  Cubic spline interpolation.
   {	
	   cs2b=cs2[i]+d2*(dcs2[i]+d2*(3.0*(cs2[i+1]-cs2[i])-2.0*dcs2[i]-dcs2[i+1]+d2*(dcs2[i]+dcs2[i+1]+2.0*(cs2[i]-cs2[i+1]))));
      opac2=dotmu[i]+d2*(ddotmu[i]+d2*(3.0*(dotmu[i+1]-dotmu[i])-2.0*ddotmu[i]-ddotmu[i+1]+d2*(ddotmu[i]+ddotmu[i+1]+2.0*(dotmu[i]-dotmu[i+1]))));
      dxe=dlxedla[i]+d2*(ddlxedla[i]+d2*(3.0*(dlxedla[i+1]-dlxedla[i])-2.0*ddlxedla[i]-ddlxedla[i+1]+d2*(ddlxedla[i]+ddlxedla[i+1]+2.0*(dlxedla[i]-dlxedla[i+1]))));
		}
	
		
	san[0] = cs2b;
	san[1] = opac2;
	san[2] = dxe;
	}	


#endif