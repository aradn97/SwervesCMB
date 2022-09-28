#ifndef RECFAST_H
#define RECFAST_H

#include <math.h>
#include "initval.h"
#include "variables.h"
#include "numericx.h"

double SwitchTbEvolution;
double fHe,HO,Nnow;
double CDB,CDB_He,CB1,CB1_He1,CB1_He2,CR,CK,CK_He,CL,CL_He,CT,Bfact,H_frac,fu,Tmat;
double x_H0,x_He0,x0;
double z_eq;


double zrec[Nz0+1],xrec[Nz0+1],dxrec[Nz0+1];


void getInit(double z) {
  double rhs;
  if (z > 8000) {
    x_H0 = 1.0;
    x_He0 = 1.0;
    x0 = 1.0+2.0*fHe;
  } else if (z > 3500) {
    x_H0 = 1.0;
    x_He0 = 1.0;
    rhs = exp( 1.5 * log(CR*tcmb/(1.0+z)) - CB1_He2/(tcmb*(1.0+z)) ) / Nnow;
    rhs *= 1.0;		// ratio of g's is 1 for He++ <-> He+
    x0 = 0.5 * ( sqrt( pow(rhs-1.0-fHe,2) + 4.0*(1.0+2.0*fHe)*rhs) - (rhs-1.0-fHe) );
  } else if (z > 2000) {
    x_H0 = 1.0;
    rhs = exp( 1.5 * log(CR*tcmb/(1.0+z)) - CB1_He1/(tcmb*(1.0+z)) ) / Nnow;
    rhs *= 4.0;		//ratio of g's is 4 for He+ <-> He0
    x_He0 = 0.5*( sqrt( pow(rhs-1.0,2) + 4.0*(1.0+fHe)*rhs ) - (rhs-1.0));
    x0 = x_He0;
    x_He0 = (x0 - 1.0)/fHe;
  } else {
    rhs = exp( 1.5 * log(CR*tcmb/(1.0+z)) - CB1/(tcmb*(1.0+z)) ) / Nnow;
    x_H0 = 0.5 * (sqrt(rhs*rhs +4.0*rhs ) - rhs );
    x_He0 = 0;
    x0 = x_H0;
  }
}



void ion(int nvariable,double z, double y[], double f[]) {
  static double Lambda = 8.2245809;
  static double Lambda_He = 51.3;

  //      the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen	
  static double a_PPB = 4.309;
  static double b_PPB = -0.6166;
  static double c_PPB = 0.6703;
  static double d_PPB = 0.5300;

  // the Verner and Ferland type fitting parameters for Helium
  //  fixed to match those in the SSS papers, and now correct
  double	a_VF = pow(10.0,(-16.744));
  double	b_VF = 0.711;
  double	T_0 = pow(10.0,0.477121); // 3K
  double	T_1 = pow(10.0,5.114); // 
    
  double	x_H = y[1];
  double	x_He = y[2];
  double	x = x_H + fHe * x_He;
  double	Tmat = y[3];
  
  double n = Nnow * pow(1.0 +z,3);
  double n_He = fHe * Nnow * pow(1.0 +z,3);
  double Trad = tcmb * (1.0 +z);

  //H = cosmos.H_0();
  double Hz = HO * sqrt(pow((1.0+z),4)/(1.0+z_eq)*(omegab+omegac) + (omegab+omegac)*pow((1.0+z),3)	+ omegav) ; // general H(z) 
  
  // Get the radiative rates using PPQ fit (identical to Hummer's table)
  double Rdown=1e-19*a_PPB*pow(Tmat/1e4,b_PPB) /(1.0+c_PPB*pow(Tmat/1e4,d_PPB));
  double Rup = Rdown * pow(CR*Tmat,1.5)*exp(-CDB/Tmat);
  
  // calculate He using a fit to a Verner & Ferland type formula
  double sq_0 = sqrt(Tmat/T_0);
  double sq_1 = sqrt(Tmat/T_1);
  //  typo here corrected by Wayne Hu and Savita Gahlaut
  double Rdown_He = a_VF/(sq_0*pow(1.0 +sq_0,1.0 -b_VF));
  Rdown_He = Rdown_He/pow(1.0+sq_1,1.0+b_VF);
  double Rup_He = Rdown_He*pow(CR*Tmat,1.5)*exp(-CDB_He/Tmat);
  Rup_He *= 4.0;	 // statistical weights factor for HeI
  //	Avoid overflow (pointed out by Jacques Roland)
  double He_Boltz;
  if (Bfact/Tmat > 680.0) He_Boltz = exp(680.0);
  else He_Boltz = exp(Bfact/Tmat);
  
  double K = CK/Hz;		// Peebles coefficient K=lambda_a^3/8piH
  double K_He = CK_He/Hz;		// Peebles coefficient for Helium
  
  //	Estimates of Thomson scattering time and Hubble time
  double timeTh=(1.0/(CT*pow(Trad,4)))*(1.0 +x+fHe)/x;	// Thomson time
  // better for hubble time :
  //double timeH=2.0/(3.0*H*pow(1.0+z,1.5)); // Hubble time
  double timeH = 2.0/(3.0*HO*pow((1.0+z),1.5)); // in s^-1

  //	calculate the derivatives
  //	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
  //	(clunky, but seems to work)
  if (x_H > 0.99) { 	// use Saha rate for Hydrogen
    f[1] = 0;
  }  else if (x_H > 0.98) {
    f[1] = (x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-CL/Tmat)) /(Hz*(1.0+z));
    // cout << "z: " << z << "  y[1] : " << y[1] << "  f[1]: " << f[1] << endl;
  } else {
    f[1] = ((x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-CL/Tmat)) *(1.0 + K*Lambda*n*(1.0-x_H)));
    f[1] /=  Hz*(1.0+z)*(1.0/fu+K*Lambda*n*(1.0-x_H)/fu +K*Rup*n*(1.0-x_H));  // was (1-x) typo thank to Jens Chluba
    // compare this typo also to correct version in RecfastAlpha
  }
  //   turn off the He once it is small
  if (x_He < 1e-15) f[2]=0; 
  else {
    f[2] = ((x*x_He*n*Rdown_He - Rup_He*(1.0-x_He)*exp(-CL_He/Tmat)) *(1.0+ K_He*Lambda_He*n_He*(1.0-x_He)*He_Boltz));
    f[2] /= Hz*(1.0+z) * (1.0 + K_He*(Lambda_He+Rup_He)*n_He*(1.0-x_He) *He_Boltz); 
  }
  //    follow the matter temperature once it has a chance of diverging
  // M.Doran:: I'm not happy with this, as it still bears the possibility
  // of switching over a little too late or risking oscillations in T_b when 
  // switching on the coupling, such that it can run to it's "true" value
  // However, we drop the first 4-5 z's after the switch and therefore
  // we will most propably not see many oscillations. In addition,
  // at the time this is happening, the baryon sound speed is already
  // negligible and we don't have to worry for the soundspeed.
  // Yet,  I still would prefer a "perfect"  T_matter evolution...
  if (timeTh < H_frac*timeH) 
    f[3]=Tmat/(1.0+z); // Tmat follows Trad
  else {
    if (SwitchTbEvolution == 0) SwitchTbEvolution = z;
    f[3]= CT * pow(Trad,4) * x / (1.0+x+fHe) * (Tmat-Trad) / (Hz*(1.0 +z)) + 2.0*Tmat/(1.0+z);
  }            

  //  cout << "curioisity: " << z << "  :: " <<  (1.0 + K*Lambda*n*(1.0-x_H)) / (1.0/fu+K*Lambda*n*(1.0-x)/fu +K*Rup*n*(1.0-x)) << endl;


}          



void recfast() 
{  //  ofstream helium("helium_old.dat");
	double zinitial = 1e4;
   double z = zinitial;
   double zfinal = 0;

	int ind,nw;				//dvark paramater
	double cw[25];	
	double w[10][nvar0+1];
   SwitchTbEvolution = 0;  // this will be the redshift at which ion() switches from tight for T_b to non-tight


   double mu_H = 1.0/(1.0-yhe);     //Mass per H atom
//  OmegaB = cosmos.omega_b();
//  OmegaT = cosmos.omega_m(); // total matter 
//  Tnow = cosmos.T_cmb();
//  double yhe = cosmos.Y_he();


//  mu_T = 4.0/(4.0-3.0*yhe); //Mass per atom
   fHe = yhe/(4.0*(1.0-yhe)); // n_He_tot / n_H_tot
 
   HO = h0/(3.0856775807e19);  // in seconds (that's originals recfast's HO)
  
  //  today's number density of hydrogen
   Nnow = 3*HO*HO*omegab/(8.0*3.141592653589*G*mu_H*m_H);
  // UNUSED double n = Nnow * pow(1 +z,3);
   double fnu = (annur*7.0/8.0)*pow(4.0/11.0, 4.0/3.0);
  //fnu = (annur*7.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
   z_eq = (3.0*(HO*C)*(HO*C)/(8.0*3.141592653589*G*a*(1.0+fnu)*pow(tcmb,4)))*(omegab+omegac);
   z_eq = z_eq - 1.0;
	
  // set up some constants so they don't have to be calculated later
  
   double Lalpha = 1.0 / L_H_alpha;
   double Lalpha_He = 1.0 / L_He_2p;
   double DeltaB = h_P*C*(L_H_ion-L_H_alpha);
   CDB = DeltaB/k_B;
   double DeltaB_He = h_P*C*(L_He1_ion-L_He_2s); // 2s, not 2p
   CDB_He = DeltaB_He/k_B;
   CB1 = h_P*C*L_H_ion/k_B;
   CB1_He1 = h_P*C*L_He1_ion/k_B; // ionization for HeI
   CB1_He2 = h_P*C*L_He2_ion/k_B; // ionization for HeII
   CR = 2.0*3.141592653589*(m_e/h_P)*(k_B/h_P);
   CK = pow(Lalpha,3)/(8.0*3.141592653589);
   CK_He = pow(Lalpha_He,3)/(8.0*3.141592653589);
   CL = C*h_P/(k_B*Lalpha);
   CL_He = C*h_P/(k_B/L_He_2s); // comes from det.bal. of 2s-1s
   CT = (8.0/3.0)*(sigma/(m_e*C))*a;
   Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B;

   double y[4]; // 4 dim for 3 variables for odeint(); 

  //    Matter departs from radiation when t(Th) > H_frac * t(H)
  //    choose some safely small number
   H_frac = 1e-3;
      
  //    Fudge factor to approximate for low z out of equilibrium effect
   fu=1.14;
      
  //   Set initial matter temperature
   y[3] = tcmb*(1.0+z);      //Initial rad. & mat. temperature
   Tmat = y[3];
  
   getInit(z); // set x_H0, x_He0 and x0
   y[1] = x_H0;
   y[2] = x_He0;
    
  //     OK that's the initial conditions, now start writing output file
   int Nz=Nz0;//10000;
   double zstart,zend;

//	Set up work-space stuff for DVERK
	ind  = 1;
	nw   = 3;
	for(int i = 0;i<=24;i++)
		cw[i] = 0.0;

//   double hnext = 1e-1; // initial stepsize guess for odeint
   for (int i = 1; i <= Nz;i++) 
   {
		//     calculate the start and end redshift for the interval at each z
      //     or just at each z
      zstart = zinitial + (i-1)*(zfinal-zinitial)/Nz;
      zend   = zinitial + i*(zfinal-zinitial)/Nz;
    
		    //     Use Saha to get x_e, using the equation for x_e for ionized helium
      //     and for neutral helium.
      //     Everything ionized above z=8000.  First ionization over by z=5000.
      //     Assume He all singly ionized down to z=3500, then use He Saha until
      //     He is 99% singly ionized, and *then* switch to joint H/He recombination.
 	   //    double dy[4];  // for getting the accurate d (ln T_baryon / dz) used only at late times
 	   //    double dTbTb = 1.0 / (1.0 + z); //   d (ln T_baryon) / dz 
      double rhs;
      z = zend;
      if (zend > 8000) 
      {
			x_H0 = 1;
      	x_He0 = 1;
      	x0 = 1.0 + 2.0*fHe;
      	y[1] = x_H0;
      	y[2] = x_He0;
      	y[3] = tcmb*(1.0+z);
  	  		} 
		else if (z > 5000) 
		{            
      	x_H0 = 1;
      	x_He0 = 1;
      	rhs = exp( 1.5 * log(CR*tcmb/(1.0+z)) - CB1_He2/(tcmb*(1.0+z)) ) / Nnow;
      	rhs *= 1.0;       // ratio of g's is 1 for He++ <-> He+
      	x0 = 0.5 * ( sqrt( pow( rhs-1.0 -fHe,2.0)  + 4.0*(1.0 +2.0*fHe)*rhs) - (rhs-1.0-fHe) );
      	y[1] = x_H0;
      	y[2] = x_He0;
      	y[3] = tcmb*(1.0+z);
    		} 
		else if (z > 3500) 
		{
      	x_H0 = 1.0;
      	x_He0 = 1.0;
      	x0 = x_H0 + fHe*x_He0;
      	y[1] = x_H0;
      	y[2] = x_He0;
      	y[3] = tcmb*(1.0 +z);
    		} 
		else if (y[2] > 0.99) 
		{
      	x_H0 = 1.0;
      	rhs = exp(1.5 * log(CR*tcmb/(1.0+z)) - CB1_He1/(tcmb*(1.0 +z)) ) / Nnow;
      	rhs *= 4.0;      //ratio of g's is 4 for He+ <-> He0
      	x_He0 = 0.5 * (sqrt( pow(rhs-1.0,2.0) + 4.0*(1.0+fHe)*rhs ) - (rhs-1.0));
      	x0 = x_He0;
      	x_He0 = (x0 - 1.0)/fHe;
      	y[1] = x_H0;
      	y[2] = x_He0;
      	y[3] = tcmb*(1.0+z);
      	//      helium << z << " " << x_He0 << endl;
    		} 
		else if (y[1] > 0.99) 
		{
      	rhs = exp( 1.5 * log(CR*tcmb/(1.0+z)) - CB1/(tcmb*(1.0+z)) ) / Nnow;
      	x_H0 = 0.5   * (sqrt( rhs*rhs+4.0*rhs ) - rhs );
      	dverk(3,nw,ion,zstart,y,zend,1.0e-5,ind,cw,nw,w);
      
			//      hnext = Miscmath::odeint(y,3,zstart,zend,1e-9,hnext,0, (moDerivs)&Recfast::ion,*this,false);  
      	y[1] = x_H0;
      	x0 = y[1] + fHe*y[2];
			//      ion(zend,y,dy); // to get derivative of Tb
			//     dTbTb = dy[3] / y[3]; 
      	//      helium << z << " " << y[2] << endl;
			} 
		else 
		{
	   	dverk(3,nw,ion,zstart,y,zend,1.0e-5,ind,cw,nw,w); 
			//      hnext = Miscmath::odeint(y,3,zstart,zend,1e-9,hnext,0, (moDerivs)&Recfast::ion,*this,false,Miscmath::bstoer);
      	x0 = y[1] + fHe*y[2];
			//      ion(zend,y,dy); // to get derivative of Tb
			//      if (SwitchTbEvolution > 0) {  // so we have switched already 
			//	if (SwitchTbEvolution - zend > 5) { // drop the first z's
			//	  dTbTb = dy[3] / y[3]; 
			} // i.e. do nothing if it just switched, keep canonical dTbTb value
//		Trad = tcmb*(1.0+zend);
	  	Tmat = y[3];
//	  	x_H = y[1];
//	  	x_He = y[2];
//	  	x = x0;

		zrec[i]=zend;
      xrec[i]=x0;		
      printf("\nzend = %e x0 = %e",zend,x0);
      }
     
     spline(zrec,xrec,Nz,1.0e40,1.0e40,dxrec); 
//   }
//    noteX->set(-zend, x0);  // set the spline 
//    noteTb->set(y[3]); // y[3] holds baryon temperature
//    noteTbDeriv->set( -(1+zend) * dTbTb); //  d ln Tb / d ln a
    // eqn (68) of Ma & Bertschinger, the + instead of minus comes from 
    // d ln Tb_ dln a = -(1+z) d ln Tb / dz 
//    double cs2 = y[3] * BoltzmannOverWeightAndC2() * (1.0 +  1.0/3.0 *(1+zend)*dTbTb);
//    noteCs2->set(cs2);
//  }  // end the loop over the z's
//  noteX->arm(Spline::all);
  //  noteX->dump("recombin");
	}


#endif