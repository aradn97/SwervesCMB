#ifndef	_VARIABLES_H_
#define	_VARIABLES_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "initval.h"

int lmo;																  // Maximum l values to be calculated for the power spectrum + 300

double omegab,omegac,omegav,omegan,omegak,omega;		  // Omega = 1-OmegaK
double omegavdyn;										           // Omega dark energy at any redshift (OmegaV * a^(3(1+w)))
double h0,tcmb,yhe,annur,annunr; 						     // Annur - Number of Neutrino specises (3.04)
                													  // Anunr - Number of Massive Neutrino specises 


double grhom,grhog,grhor,grhonr;								  // g -> 8.pi.G/c^2		
// grhom  - (8*pi*G/c^2)*rhocr = 3*Ho^2/c^2
// grhog  - (8*pi*G/c^2)*(4*Sigma_sb/c^3)*T^4           // for photons 
// grhor  - (7/8)*((4/11)^(4/3))*grhog                  // for massless neutrinos
// grhonr - (7/8)*((43/(11*g))^(4/3))*grhog             // for massive neutrinos


double akmax0;															//	Maximum k*tau0 for integration. User input for "ketamax".			
double akmaxt;															// Maximum k till which transfer functions are required			
double armax;															// a where visibility function is maximum		
double adotrad;														// da/dtao in the radiation dominated era

double tau0;															// Conformal time now
double taurst;															// Recombination start time
double taurmax;														// Where visibility function is maximum
double taur;															// Conformal time of the recombination era (maximum visibility)	

double dtaurec;														// Time step during recombination ... it is modified from cmbfalst() and finithermo() 

int nri0,nri; 														   // nri0 <- Reionization 1st step number of discritization
																			// nri <- Reionization 2nd step number of discritization
int j2ri1;																// j2ri1 <- timestep at (9/10)taurist where taurist the time when reionization had started 
int jrmax;														   	// point where visibility function is maximum	


int lmaxg,lmaxnr;														// lmaxg <- Photon number of multipoles
																			// lmaxnr <- Nutrino (massless) number of multipoles
int iq0,iq1,iq2;														// Used for fixing number of equations
int lmaxnu,nqmax;														//	Num of discritization of massive neutrino multipoles and number of momentum step	


double rhonu,drhonu;													// rhonu <- Massive Neutrino density
																			// drhonu <- Derivative of Massive Neutrino density
double dlfdlq[nqmax0+1];											// dlfdlq <- Massive neutrino d_log_f/d_log_q where f is the neutrino distribution function
double qdn[nqmax0+1];												// q^3 dq / (1+exp(q)) ..... used by nu2()

double epsw;															// This will set a limit on webnumber

double adotoa;															// a_dot/a ......... it is used by different functions 

double hdot,dgshear;													// hdot <- scalar metric perturbation
																			// dgshear <- derivative of share of photons			

double aksplit;														// Wave number corresponding the horizon scale

double dalphansdlnk[nnmax+1]={0};										// 
double an[nnmax+1],alphans[nnmax+1];							//	Scaler spectral index														
																									
																		
/****************************************************************************************/
//		Flag Variables																							 //
/****************************************************************************************/

int dimflag;											           // Flag variable  1 - 5th dim,  
																	 	  //                0 - 4 dim

int rcflag;															  // Flag variable 0 - peebles
																		  // 					 1 - recfast																	 	  

int lensflag; 														  // Flag lensing/not  0 - for unlensed Cl 
																		  //                   1 - for lensed Cl linear evalution
																		  //					     2 - for lensed Cl nonlinear evalution

int initfl;  														  //	Flag initial condition  1 - Adiabatic Initial Condition	 
																		  // 									2 - Isothermal CDM
																		  //									3 - Isothermal Baryonic
																		  //									4 - Isothermal seed Condition

int kcutflag;														  //	K Split			 1 - Only k > horizon scale are requested 
																		  // 						-1	- Only k < horizon scale are requested		

int ilin;															  // Flag variable for lensing : 0 - Linear Evolution
																		  // 										1 - Nonlinear Evolution		
																		  		
unsigned int twofield;												// Flag variable for inflationary model :	0 - Single field standard inflation
																			//													 :	1 - Two field model

unsigned int powerfileflag;										// Flag variable for twofield power model : 0 - Use standard function of this code
																			//														: 1 - Read powers from a file 
				 
/****************************************************************************************/
//		Storing Power spectrums																				 //
/****************************************************************************************/

double cl[lmax+1][nnmax+1]={{0}};							// 
double clpr[lmax+1][nnmax+1]={{0}};							// 
double cpl[lmax+1][nnmax+1]={{0}};							// Scaler  Cl's		
double cplpr[lmax+1][nnmax+1]={{0}};						// 
double ccl[lmax+1][nnmax+1]={{0}};							//	 
double cclpr[lmax+1][nnmax+1]={{0}};						// 
double ckkl[lmax+1][nnmax+1]={{0}};							// 
double ckklpr[lmax+1][nnmax+1]={{0}};						// 
double ctkl[lmax+1][nnmax+1]={{0}};							// 
double ctklpr[lmax+1][nnmax+1]={{0}};						//	 

double cladb[lmax+1][nnmax+1]={{0}};						// 
double clpradb[lmax+1][nnmax+1]={{0}};						// 
double cpladb[lmax+1][nnmax+1]={{0}};						// Scaler  Cl's	- adiabatic two field model	
double cplpradb[lmax+1][nnmax+1]={{0}};					// 
double ccladb[lmax+1][nnmax+1]={{0}};						//	 
double cclpradb[lmax+1][nnmax+1]={{0}};					// 

double cliso[lmax+1][nnmax+1]={{0}};						// 
double clpriso[lmax+1][nnmax+1]={{0}};						// 
double cpliso[lmax+1][nnmax+1]={{0}};						// Scaler  Cl's	- isocurvature two field model	
double cplpriso[lmax+1][nnmax+1]={{0}};					// 
double ccliso[lmax+1][nnmax+1]={{0}};						//	 
double cclpriso[lmax+1][nnmax+1]={{0}};					// 

double clcross[lmax+1][nnmax+1]={{0}};						// 
double clprcross[lmax+1][nnmax+1]={{0}};					// 
double cplcross[lmax+1][nnmax+1]={{0}};					// Scaler  Cl's	- Cross two field model	
double cplprcross[lmax+1][nnmax+1]={{0}};					// 
double cclcross[lmax+1][nnmax+1]={{0}};					//	 
double cclprcross[lmax+1][nnmax+1]={{0}};				 	// 

double cltotal[lmax+1][nnmax+1]={{0}};						// 
double clprtotal[lmax+1][nnmax+1]={{0}};					// 
double cpltotal[lmax+1][nnmax+1]={{0}};					// Scaler  Cl's	- Total two field model	
double cplprtotal[lmax+1][nnmax+1]={{0}};				 	// 
double ccltotal[lmax+1][nnmax+1]={{0}};					//	 
double cclprtotal[lmax+1][nnmax+1]={{0}};				 	// 

/****************************************************************************************/
//		Storing Power spectrums	(Tensor)																	 //
/****************************************************************************************/
double cltt[l0max+1][nnmax+1]={{0}};						// T mode
double clet[l0max+1][nnmax+1]={{0}};						// E mode
double clbt[l0max+1][nnmax+1]={{0}};						// B mode
double clct[l0max+1][nnmax+1]={{0}};						// TE mode

/****************************************************************************************/
//		Reionization Variables																				 //
//    Developer's Guide : (This part is rellated to reion.h file) 				          //
/****************************************************************************************/
double zri,taurist,zristp,tauristp;																		 	
double rif,optdlss;																							 	

// zri      - Reionization redshift																		 //	
// zristp   - Redshift at which reionization stop (Just for making grid)					 //		
// taurist  - Conformal time at reionization															 //	
// tauristp - Conformal time at which reionization stop
// rif      - reionization fraction																		 //	
// optdlss  - optical depth to lss 																		 //	

//////////////////////////////////////////////////////////////////////////////////////////



/****************************************************************************************/
//   These variables are calculated from the recombination process							 //
/****************************************************************************************/	

// These variables are from the initial greading
double cs2[nthermo+1];								// Baryon sound speed square 
double dotmu[nthermo+1];							// dmu/deta ... scattering rate
double sdotmu[nthermo+1] = {0};					// Int ((dmu/deta)deta,0,ith time)  .... mu i.e. integration of dmu/deta
double dlxedla[nthermo+1];							// d ln(Xe) / d ln(a) 	
double emmu[nthermo+1];								// exp(-mu)
double ddotmu[nthermo+1];							// Derivative and double derivative of the abive quantities	
double dddotmu[nthermo+1];							//             "   "
double ddddotmu[nthermo+1];						//             "   "
double demmu[nthermo+1];							//					"   "
double ddlxedla[nthermo+1];						//					"   "			  
double dcs2[nthermo+1];								// 				" 	 "	 


// These variables are after new grid is defined. These are the values of the same variables as the above list but interpolated in the new time grid
double vis[nstep0+1],dvis[nstep0+1];			// vis <- visibility function	g = (dmu/dtau)*exp(-mu)									
double ddvis[nstep0+1],expmmu[nstep0+1];		// opac <- 	 dmu/deta ... scattering rate ...... interpolated	 
double opac[nstep0+1],dopac[nstep0+1];			// expmmu <- exp(-mu) ....... interpolated emmu
															// Others are their derivatives and 2nd derivatives																		
double dtau1[nstep0+1],dtau2[nstep0+1];		// dtau1 <- Tau(i+1)-Tau(i)  ............     time step at each grid point 
															// dtau2 <- (Tau(i+1)-Tau(i-1))/2  ......     time step at each grid point
double atau0[nstep0+1]={0};						// atau0 <- time steps in the new grid 


/****************************************************************************************/
//   These variables are calculated from the recombination process							 //
/****************************************************************************************/	
double ak,ak2;											//
double apowers;

double clts[l0max+1][nnmax+1],cles[l0max+1][nnmax+1];
double clkk[l0max+1][nnmax+1],cltk[l0max+1][nnmax+1];
double clcs[l0max+1][nnmax+1],clbs[l0max+1][nnmax+1];

double cltsadb[l0max+1][nnmax+1],clesadb[l0max+1][nnmax+1],clcsadb[l0max+1][nnmax+1];
double cltsiso[l0max+1][nnmax+1],clesiso[l0max+1][nnmax+1],clcsiso[l0max+1][nnmax+1];
double cltscross[l0max+1][nnmax+1],clescross[l0max+1][nnmax+1],clcscross[l0max+1][nnmax+1];
double cltstotal[l0max+1][nnmax+1],clestotal[l0max+1][nnmax+1],clcstotal[l0max+1][nnmax+1];

int nn;													
int ntf;

double anorm[ntfmax+1][nnmax+1];
	

#endif