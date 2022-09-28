#ifndef	_INITVAL_H_
#define	_INITVAL_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


// NEOTRINOS (massive)

// For initiallising the massive neutinos. 
// For details check (P3)

#define constan 5.68219698     // 7(Pi^4)/120
#define zeta3 1.20205690       // Zeta(3) for calculating the number densities of fermions



// NEUTRINOS (massive)

// For details check (P3)

#define lm3 4                  // Number of discritization of k           (Neutrinos massive)
#define nql 15                 // Number of discritization of the energy  (Meutrinos massive)

//For details check initnul()

#define nqmax0 15					 // Number of discritization of the integration x^3 dx /(1+exp(x))
#define nrhopn 10000				 // Number of discritization of  pressure and the density (Neutrinos massive) 

// used by finithermo() in recombination.h
#define nthermo 10000			 // Number of points at which opacity and the baryon temperatures are calculated
#define nstep0 7000				 // Number of discritizations of the time step (after recombination for the actual calculation)


#define lmax0 12					 //Number of photon equations for each polarizations
#define lmaxnr0 25				 //Number of relativistic neutrino equations 
#define nqmax0 15					 //Number of q's at which massive neutrino equations will be written
#define lmaxnu0 25				 //Number of massive neutrino equations for each q. 	


#define l0max 5300             //l0max should be more than the maximum l by atleast 300
#define nnmax 11					 //Maximum number of spectral indexes can be calculated in one run. Change if needed
#define ntfmax 211				 //Maximum number of transfer functions can be calculated in one run. Change if needed

#define nvar0 9+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1)+1	// Total number of equations

#define Nz0 10000

#define lmax 20+l0max/10       //Number of Cl's to be calculated. lmax slould be more than the dimension of the l[] array. 
/***********************************************************************************/
// Constant Values 
// These values are used by recfast function
/***********************************************************************************/

//	C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
//	m_e,m_H: electron mass and mass of H atom in SI

//	sigma: Thomson cross-section
//	a: radiation constant for u=aT^4


//	L_H_ion: level for H ionization in m^-1
//	L_H_alpha: level for H Ly alpha in m^-1
//	L_He1_ion: level for HeI ionization
//	L_He2_ion: level for HeII ionization
//	L_He_2s: level for HeI 2s
//	L_He_2p: level for HeI 2p
//	Lalpha: Ly alpha wavelength in SI
//	Lalpha_He: Helium I 2p-1s wavelength in SI

double C  = 2.99792458e8;
double k_B = 1.380658e-23;
double h_P = 6.6260755e-34;
double m_e = 9.1093897e-31;
double m_H = 1.673725e-27;

double sigma = 6.6524616e-29;
double a = 7.565914e-16;
  
double G = 6.67259e-11;

//   2 photon rates and atomic levels in SI units

double L_H_ion = 10967877.37;			//level for H ion. (in m^-1)
double L_H_alpha = 8225916.453;		//averaged over 2 levels
double L_He1_ion = 19831077.2;		//from Drake (1993)
double L_He2_ion = 43890888.63;		//from JPhysChemRefData (1987)
double L_He_2s = 16627743.4;			//from Drake (1993)
double L_He_2p = 17113489.1;			//from Drake (1993)

#endif