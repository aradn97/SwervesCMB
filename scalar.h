#ifndef SCALAR_H_
#define SCALAR_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "initval.h"
#include "variables.h"
#include "numericx.h"
#include "darkenergy.h"
#include "neutrino.h"
#include "others.h"
#include "reion.h"
#include "recombination.h"
#include "curvature.h"

#include "VPC.h"

class scalar
{

   double adotoa;        // a_dot/a ......... it is used by different functions
   double hdot, dgshear; // hdot <- scalar metric perturbation
                         // dgshear <- derivative of share of photons

   double rhonu, drhonu; // rhonu <- Massive Neutrino density
                         // drhonu <- Derivative of Massive Neutrino density

 public:
   /************************************************************************************************************************************/
   /*		It will store the transfer function at different point																								*/
   /*		It will be used by the function cmbflat()																													*/
   /************************************************************************************************************************************/

   void outtransf(int n, double y[], double curv, int itf)
   {
      double tfc, tfb, tfg, tfn, a, deltan, tfnm, beta, pnu, rhonu, xsanpass[2];
      double zsanpass[4];
      double drhonu;

      neutrino Neutrino;
      beta = Variables.ak; //sqrt(ak-curv);
      Variables.ak2 = beta * beta;

      tfc = y[4] / Variables.ak2;
      tfb = y[6] / Variables.ak2;
      tfg = y[8] / Variables.ak2;
      tfn = y[10 + 2 * Variables.lmaxg] / Variables.ak2;


      
      if (Variables.amnu != 0.0)
      {
         a = y[1];
         Neutrino.nu1(a, xsanpass);
         rhonu = xsanpass[0];
         pnu = xsanpass[1];
         Neutrino.nu2(a, zsanpass, y + Variables.iq0, y + Variables.iq1, y + Variables.iq2);
         drhonu = zsanpass[0];

         deltan = drhonu / rhonu;
         tfnm = deltan / Variables.ak2;
      }

      GlobalArray.power[6 * GlobalArray.istore[itf] + 1][itf] = beta / (Variables.h0 / 100.0);
      GlobalArray.power[6 * GlobalArray.istore[itf] + 2][itf] = tfc;
      GlobalArray.power[6 * GlobalArray.istore[itf] + 3][itf] = tfb;
      GlobalArray.power[6 * GlobalArray.istore[itf] + 4][itf] = tfg;
      GlobalArray.power[6 * GlobalArray.istore[itf] + 5][itf] = tfn;

      if (Variables.amnu == 0.0)
         GlobalArray.power[6 * GlobalArray.istore[itf] + 6][itf] = 0.0;
      else
         GlobalArray.power[6 * GlobalArray.istore[itf] + 6][itf] = tfnm;

      GlobalArray.istore[itf] = GlobalArray.istore[itf] + 1;
      
      //printf("\n111 %e %e %e %e %e %e",beta/(Variables.h0 / 100.0),Variables.ak2,tfc,tfb,tfg,tfn);      
   }

   /************************************************************************************************************************************/
   /* This part is based on M.Zaldarriaga's Ph.D. Thesis																											*/
   /* It will calculate the source term for the CMBR 																												*/
   /************************************************************************************************************************************/

   /************************************************************************************************************************************/
   /*		finital will calculate the initial values for the source terms depending on the choosen initial condition							*/
   /*		It will use the functions : dynrho(),																														*/
   /************************************************************************************************************************************/

   void finitial(double y[], double tau)
   {
      double rhosan[2];
      double rhonu, pnu;
      double a, a2;

      double grho, fracnu, s, gpres;
      double psi, C, akt2, h, deltac, deltag, deltab, deltar;
      double thetac, thetag, thetab, thetar, shearr, ahdot;
      double dq, q, aq, v, akv, f1, yrad, weos;
      double rc_phi, rc_psi, expq, eta;
      double deltan, thetan, delta0;
      double dlfdlq1;

      int ind;

      double ak, ak2;
      double beta, beta2;

      VPC vpc;

      ak = Variables.ak;
      ak2= Variables.ak2;

      neutrino Neutrino;
      darkenergy Darkenergy;
      curvature Curvature(-Variables.omegak);   // Need to change XYZ

	   Variables.aux = 1;      
      
      if(fabs(Variables.omegak)>0.001){
         //curv=-omegak*h0*h0/csp/csp;
         beta = ak;
         beta2 = beta*beta;
   	   ak=sqrt(beta2-Variables.curv);
	      ak2=ak*ak;

	      Variables.aux=1.0/(1.0-3.0*Variables.curv/ak2);
      }

      //printf("\n %e %e",tau, Variables.adotrad);
      //exit(1);

      a = tau * Variables.adotrad; // a at the initial point
      a2 = a * a;


      Neutrino.nu1(a, rhosan); // pressure and density at initial point
      rhonu = rhosan[0];
      pnu = rhosan[1];

      //printf("%e %e %e %e",a,rhonu,pnu,ak2);
      //exit(1);

      Variables.omegavdyn = Variables.omegav * Darkenergy.dynrho(a);

      if(Variables.VPCflag == 1)  
         Variables.omegavdyn = Variables.omegavdyn * vpc.Variable_Lambda(a);


      weos = Darkenergy.wdyn_func(a);

      if (Variables.dimflag == 0) // 8*pi*G*rho*a**2. energy density for FRW and 5 dimensions
         grho = vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + vpc.Variable_grhonr(a) * Variables.annunr * rhonu) / a2 + vpc.Variable_grhom(a) * Variables.omegavdyn * a2;
      else
      {
         grho = sqrt(grho) + sqrt(Variables.omegav * vpc.Variable_grhom(a) * a2);
         grho = grho * grho;
      }
      // pressure 8*pi*G*P*a**2.
      gpres = ((vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur) / 3.0 + vpc.Variable_grhonr(a) * Variables.annunr * pnu) / a2 + weos * vpc.Variable_grhom(a) * Variables.omegavdyn * a2;

      s = grho + gpres;
      fracnu = (4.0 / 3.0) * (vpc.Variable_grhor(a) * Variables.annur + vpc.Variable_grhonr(a) * Variables.annunr) / (a2 * s); // Nutrino Fraction

      // Use yrad=rho_matter/rho_rad to correct initial conditions for matter+radiation.
      yrad = vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) * a / (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + vpc.Variable_grhonr(a) * Variables.annunr * rhonu);

      //printf("\nHHHHHHHHHHHHHHHHHHHHHHH\n %e %e %e %e %e %e",a, grho,gpres,s,fracnu,yrad);
      //exit(1);

      /**********************************************************************************************************/
      /* Choose one of the following four cases for initial conditions, or add your own.                        */
      /**********************************************************************************************************/

      /***********************************************************************************************************/
      /*  First case.												   */
      /*  Isentropic ("adiabatic") initial conditions. normalize to zeta=1					   */
      /*  Check the equation set (96),(97),(98) Ma,Bertschinger                           			   */
      /***********************************************************************************************************/
      
      if (Variables.initfl == 1)
      {
         if(fabs(Variables.omegak)<0.001)
            psi = -1.0 / (1.5 + 2.0 * fracnu / 5.0); //fracnu = rhonu/(rhogama+rhonu)	Eq(98f)
         else
            psi = -1.0;
         C = (15.0 + 4.0 * fracnu) / 20.0 * psi;
         //printf("\n %e %e %e",C,fracnu,psi);exit(1);
         akt2 = (ak * tau) * (ak * tau);
         h = C * akt2 * (1.0 - 0.2 * yrad);                                                                    //Eq(96g)
         eta = 2.0 * C - (5.0 + 4.0 * fracnu) * C * akt2 * (1.0 - yrad / 3.0) / (6.0 * (15.0 + 4.0 * fracnu)); //Eq(96h)
         f1 = (23.0 + 4.0 * fracnu) / (15.0 + 4.0 * fracnu);                                                   //Eq(96e f1=thetanu/thetagamma)
         deltac = -0.5 * h;                                                                                    //Eq(96a, 96b, 96g)
         deltag = -2.0 * h * (1.0 - akt2 / 36.0) / 3.0;                                                        //Eq(96a)
         deltab = 0.75 * deltag;                                                                               //Eq(96b)
         deltar = -2.0 * h * (1.0 - akt2 * f1 / 36.0) / 3.0;                                                   //Eq(96a, 96b)
         thetac = 0.0;                                                                                         //Eq(96c)
         thetag = -C * akt2 * akt2 / (tau * 18.0);                                                             //Eq(96d)
         thetab = thetag;                                                                                      //Eq(96d)
         thetar = f1 * thetag;                                                                                 //Eq(96e)
         shearr = (4.0 / 15.0) * (ak2 / s) * psi * (1.0 + 7.0 / 36.0 * yrad)/Curvature.b(2,ak2);                        //Eq(98d)
         ahdot = 2.0 * C * ak2 * tau * a * (1.0 - 0.3 * yrad);                                       //Eq(94a) a*hdot
         //printf("\n %e",Variables.omegak);
         //printf("\n %e %e %e %e %e %e %e %e %e %e",psi,C,akt2,h,eta,f1,deltac,deltag,deltab,deltar);
         //printf("\n %e %e %e %e %e %e",thetac,thetag,thetab,thetar,shearr,ahdot);
         //printf("\n XXXXXXXXXXXXXXXXXXXXXX");
         //exit(1);
      }

      /***********************************************************************************************************/
      /*  Second case.																														  */
      /*  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.												  */
      /*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois	
      /*     and "The General Primordial Cosmic Perturbation" - M. Bucher
      /*     and "Isocurvature initial conditions for second order Boltzmann solvers" by P. Carrilho									  */
      /*  Check Eq (7),(8),(9)     Sb=0    Snu=0																					  */
      /***********************************************************************************************************/
      else if (Variables.initfl == 2)
      {
         delta0 = 1.0;
         h = delta0 * yrad * (1.0 / (1.0 + Variables.omegab / Variables.omegac) - 0.5 * yrad);
         //printf("\n%e %e",yrad,h);
         //exit(1);
         deltac = delta0 - 0.5 * h; // Perturb CDM
         deltag = -2.0 / 3.0 * h;   //
         deltab = 0.75 * deltag;    // Sb = 0
         deltar = deltag;           // Snu = 0
         thetac = 0.0;
         thetag = -h * ak2 * tau / 12.0;
         thetab = thetag; //Eq(11, 7)
         thetar = thetag; //Eq(11)
         shearr = 0.0;
         ahdot = Variables.adotrad * h * (1.0 - 0.5 * yrad);
         eta = -h / 6.0;
         if(fabs(Variables.omegak)>0.001)
            eta=(eta+0.5*Variables.curv*h/ak2);
         //printf("\n  %e %e %e %e %e %e",h,eta,deltac,deltag,deltab,deltar);
         //printf("\n %e %e %e %e %e %e",thetac,thetag,thetab,thetar,shearr,ahdot);
         //exit(1);            
      }

      /***********************************************************************************************************/
      /*  Third case.																														  */
      /*  Isocurvature Baryon initial conditions: perturb only baryons as a --> 0.											 */
      /*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois	
      /*     and "The General Primordial Cosmic Perturbation" - M. Bucher
      /*     and "Isocurvature initial conditions for second order Boltzmann solvers" by P. Carrilho								  */
      /*  Check Eq (7),(8),(9)     Sr=0    Snu=0					   															  */
      /***********************************************************************************************************/
      else if (Variables.initfl == 3)
      {
         delta0 = 1.0;
         h = delta0 * yrad * (1.0 / (1.0 + Variables.omegac / Variables.omegab) - 0.5 * yrad);
         deltab = delta0 - 0.5 * h;
         deltac = -0.5 * h;
         deltag = -2.0 / 3.0 * h;
         deltar = deltag;
         thetac = 0.0;
         thetag = -h / 12.0 * ak2 * tau;
         thetab = thetag;
         thetar = thetag;
         shearr = 0.0;
         ahdot = Variables.adotrad * h * (1.0 - 0.5 * yrad);
         eta = -h / 6.0;
         if(fabs(Variables.omegak)>0.001)
            eta=(eta+0.5*Variables.curv*h/ak2);

         //printf("\n %e %e %e %e %e %e",h,eta,deltac,deltag,deltab,deltar);
         //printf("\n %e %e %e %e %e %e",thetac,thetag,thetab,thetar,shearr,ahdot);
         //exit(1);                                 
      }
      /***********************************************************************************************************/
      /*  Fourth case.																														  */
      /*  Isocurvature Seed initial conditions: unperturbed everything 													     */
      /*	 Check "Isocurvature cosmological perturbations and the CMB" - David Langlois	
      /*     and "The General Primordial Cosmic Perturbation" - M. Bucher	
      /*     and "Isocurvature initial conditions for second order Boltzmann solvers" by P. Carrilho									  */
      /*  Check Eq (7),(8),(9)     Sb=0    Snu=0	Sr =9									    									  */
      /***********************************************************************************************************/
      else if (Variables.initfl == 4)
      {
         delta0 = 1.0;
         h = delta0 * yrad * (1.0 / (1.0 + Variables.omegac / Variables.omegab) - 0.5 * yrad);
         deltab = -0.5 * h;
         deltac = -0.5 * h;
         deltag = -2.0 / 3.0 * h;
         deltar = deltag;
         thetac = 0.0;
         thetag = -h / 12.0 * ak2 * tau;
         thetab = thetag;
         thetar = thetag;
         shearr = 0.0;
         ahdot = Variables.adotrad * h * (1.0 - 0.5 * yrad);
         eta = -h / 6.0;
         if(fabs(Variables.omegak)>0.001)
            eta=(eta+0.5*Variables.curv*h/ak2);

         //printf("\n %e %e %e %e %e %e",h,eta,deltac,deltag,deltab,deltar);
         //printf("\n %e %e %e %e %e %e",thetac,thetag,thetab,thetar,shearr,ahdot);
         //exit(1);           
      }

      /***********************************************************************************************************/
      /*	 If some other value then neglect and exit form the program															  */
      /***********************************************************************************************************/
      else
      {
         printf("Its not possible. Initial condition is outside the range of values");
         exit(1);
      }

      deltan = deltar;
      thetan = thetar;

      y[1] = a;
      y[2] = ahdot;
      y[3] = eta;
      //CDM.
      y[4] = deltac;
      y[5] = thetac;
      //Baryons.
      y[6] = deltab;
      y[7] = thetab;
      //Photons (total intensity and polarization).
      y[8] = deltag;
      y[9] = thetag;

      for (int lsa = 1; lsa <= Variables.lmaxg; lsa++)
      {
         y[9 + lsa] = 0.0;                   // Ma Bertschingger Eq(64) F_gamma
         y[9 + Variables.lmaxg + lsa] = 0.0; // Ma Bertschingger Eq(64) G_gamma
      }
      // Massless neutrinos.
      y[10 + 2 * Variables.lmaxg] = deltar;
      y[11 + 2 * Variables.lmaxg] = thetar;
      y[12 + 2 * Variables.lmaxg] = shearr * 2.0; // Ma Bertschingger Eq(49) F_nu_2

      for (int l = 3; l <= Variables.lmaxnr; l++)
         y[10 + 2 * Variables.lmaxg + l] = 0.0; // Ma Bertschingger Eq(49) F_nu_l

      //  Massive neutrinos.
      dq = 1.0;
      for (int i = 1; i <= Variables.nqmax; i++)
      {
         q = i * dq - 0.5;
         aq = a * Variables.amnu / q;
         v = 1.0 / sqrt(1.0 + aq * aq);
         akv = ak * v;
         expq = exp(-q);
         dlfdlq1 = -q / (1.0 + expq);
         y[Variables.iq0 + i - 1] = -0.25 * dlfdlq1 * deltan;        // Ma Bertschingger Eq(97) Psi_nu_0 for different velovities of neutrino (massive)
         y[Variables.iq1 + i - 1] = -dlfdlq1 * thetan / (3.0 * akv); // Ma Bertschingger Eq(97) Psi_nu_1 for different velovities of neutrino (massive)
         y[Variables.iq2 + i - 1] = -0.5 * dlfdlq1 * shearr;         // Ma Bertschingger Eq(97) Psi_nu_3 for different velovities of neutrino (massive)
         for (int l = 3; l <= Variables.lmaxnu; l++)
         {
            ind = 10 + 2 * Variables.lmaxg + Variables.lmaxnr + i + l * Variables.nqmax;
            y[ind] = 0.0; // Ma Bertschingger Eq(97) (56) Psi_nu_l for diff velovities of neutrino (massive)
         }
      }
      //  Check energy constraint equation.

      if (Variables.ndyn == 1 || Variables.ndyn == 2) // Perturbations in dark energy
      {
         rc_phi = 0.0;
         rc_psi = 0.0;
         y[nvar0 - 1] = rc_phi;
         y[nvar0] = rc_psi;
      }
   }

   /************************************************************************************************************************************/
   /* 	This will calculate the derivatives of the scaler perturbation terms. This part is totally based on Ma & Bertschunger			*/
   /************************************************************************************************************************************/

   void fderivsCurved(int n, double x, double y[], double yprime[])
   { // n store total number of values in y
      double ep0 = .01;
      double xsanpass[3], ep, tau, epstc, dpnu, grho, gpres;
      double deltac, thetac, deltab, thetab, deltag, thetag, shearg;
      double polter, deltar, thetar, shearr, a2, cs2, dlxedla; //,opec
      double a, eta, rc_phi, rc_psi, photbar, pb43, ysanpass[2];
      double tcp1, tcp2, tcp3, drag, deltabdot;
      double rgrho, slip, pnu, fnu, shearnu;
      double dgrho, dgpres, dgtheta, etadot;
      double deltacdot, thetacdot, deltagdot, thetabdot;
      double rhonudot, shearnudot, adotdota, weos;
      double thetagdot, deltardot, thetardot, Fnu2dot;
      double dgprs_phi, dahdotdtau, sgrhooa2;
      double opac1, dgrho_phi, dgtheta_phi, dq, q, v, aq;
      double akv[nqmax0 + 1]; //,rc_dphi,rc_dpsi;
      double zsanpass[4];
      double tcp = 0.000000000;

      int ind, tcpa, tcpb;
      double ak,ak2,beta;

      neutrino Neutrino;
      darkenergy Darkenergy;
      numericx Numericx;
      recombination Recombination;
      curvature Curvature=curvature(-Variables.omegak);
      VPC vpc;

      ak = Variables.ak;
      ak2= Variables.ak2;

      if(fabs(Variables.omegak)>0.0001)
      {
         ak = sqrt(ak2 - Variables.curv);
         ak2 = ak*ak;
      }

      if(fabs(Variables.omegak) < 0.0001)
      {
         if (Variables.ak > Variables.epsw) // ep is used to stop the tight coupling approximation.
            ep = ep0;
         else
            ep = 0.5 * ep0;
         }
      else
      {
         if (Variables.ak > 0.6*Variables.epsw) // ep is used to stop the tight coupling approximation.
            ep = ep0;
         else
            ep = 1.1 * ep0;
         }
      

      tau = x; // time of the derivative

      a = y[1];
      eta = y[3];

      deltac = y[4]; //  CDM.
      thetac = y[5];

      deltab = y[6]; // Baryons.
      thetab = y[7];

      deltag = y[8]; //  Photons.
      thetag = y[9];
      shearg = y[10] / 2.0;

      polter = y[10] + y[9 + Variables.lmaxg] + y[11 + Variables.lmaxg]; //  Polarization term. Ma Eq(60). F2+G0+G2
      //printf("\nTe: %e %e %e %e", polter, y[10], y[9 + Variables.lmaxg], y[11 + Variables.lmaxg]);
      deltar = y[10 + 2 * Variables.lmaxg]; //  Massless neutrinos.
      thetar = y[11 + 2 * Variables.lmaxg];
      shearr = y[12 + 2 * Variables.lmaxg] / 2.0;

      a2 = a * a;
      
      //printf("\nI4 am now here ( tauminn = %e )...",Variables.tauminn);      
      Recombination.thermo(tau, xsanpass); // Get the Value of dlxedla


      cs2 = xsanpass[0];
      opac1 = xsanpass[1];
      dlxedla = xsanpass[2];



      photbar = vpc.Variable_grhog(a) / (vpc.Variable_grhom(a) * Variables.omegab * a); // Photon mass density over baryon mass density.
      pb43 = (4.0 / 3.0) * photbar;

      tcp = 0.00000000000; // Tight Coupling parameters
      tcp1 = Variables.ak / opac1;
      tcp2 = 1.0 / (opac1 * tau);
      tcp3 = ak * tau;

      //printf("\nJJJJJJJJJ111Jh %d %d %d %d %e",(tcp1 > ep),(tcp2 > ep),(tcp1 > 0.1*ep),( (tcp2 > ep) && (tcp1 > 0.1*ep),tcp ));

      if((tcp1>ep) || ((tcp2>ep) && (tcp1>(0.1*ep))))
         tcp=1.0;
      //printf("\nJJJJJJJJ222Jh %d %d %d %d %e",(tcp1 > ep),(tcp2 > ep),(tcp1 > 0.1*ep),( (tcp2 > ep) && (tcp1 > 0.1*ep),tcp ));



      tcpa = 0;
      tcpb = 0;

      epstc = Numericx.max(tcp1, tcp2);

      if(fabs(Variables.omegak)>0.0001)
      {
	      if ((tcp1 > ep) || ( (tcp2 > ep) && (tcp1 > 0.1*ep) ))
      	   tcp=1.0;
         }



      if ((epstc > 5.0e-3) && (a > 1.35e-5))
         tcpa = 1;

      epstc = epstc / (1.0 + pb43);

      if ((epstc > 5.0e-3) && (a > 1.35e-4))
         tcpb = 1;

      /*********************************************************************************************************************************/
      /*  Compute expansion rate (a dot = sqrt(8 pi G rho / 3)).																								*/
      /*********************************************************************************************************************************/

      if (Variables.amnu == 0.0) // No massive Neutrino, so no perturbation
      {
         rhonu = 1.0;
         pnu = 1.0 / 3.0;
         drhonu = 0.0;
         fnu = 0.0;
         dpnu = 0.0;
         shearnu = 0.0;
      }
      else // If massive Neutrino present then perturbe
      {
         Neutrino.nu1(a, ysanpass);
         rhonu = ysanpass[0];
         pnu = ysanpass[1];
         Neutrino.nu2(a, zsanpass, y + Variables.iq0, y + Variables.iq1, y + Variables.iq2);
         drhonu = zsanpass[0];
         fnu = zsanpass[1];
         dpnu = zsanpass[2];
         shearnu = zsanpass[3];
      }

      Variables.omegavdyn = Variables.omegav * Darkenergy.dynrho(a);

      if(Variables.VPCflag == 1)  
         Variables.omegavdyn = Variables.omegavdyn * vpc.Variable_Lambda(a);

      weos = Darkenergy.wdyn_func(a);

																							//  8*pi*G*rho*a^2 and 8*pi*G*P*a^2.

      grho = vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + Variables.annunr * rhonu * vpc.Variable_grhonr(a)) / a2 + vpc.Variable_grhom(a) * Variables.omegavdyn * a2 + vpc.Variable_grhom(a) * Variables.omegak;

      if (Variables.dimflag == 1)
         grho = pow((sqrt(grho) + sqrt(Variables.omegav * vpc.Variable_grhom(a) * a2)), 2); // energy density for  5 dimensions

      adotoa = sqrt(grho / 3.0);
      yprime[1] = adotoa * a;

      gpres = ((vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur) / 3.0 + vpc.Variable_grhonr(a) * Variables.annunr * pnu) / a2; // 8*pi*G*rho*a^2 and 8*pi*G*P*a^2.
      gpres = gpres + weos * vpc.Variable_grhom(a) * Variables.omegavdyn * a2;


      if (Variables.ndyn == 1 || Variables.ndyn == 2) // dark energy perturbations
      {
         rc_phi = y[n - 1];
         rc_psi = y[n];

         Darkenergy.dyn_nrg(a, grho, gpres, rc_phi, rc_psi, xsanpass);
         dgrho_phi = xsanpass[0];
         dgprs_phi = xsanpass[1];
         dgtheta_phi = xsanpass[2];
      }
      // Evaluate metric and massive neutrino perturbations.
      //  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
      dgrho = vpc.Variable_grhom(a) * (Variables.omegac * deltac + Variables.omegab * deltab) / a + (vpc.Variable_grhog(a) * deltag + vpc.Variable_grhor(a) * Variables.annur * deltar + vpc.Variable_grhonr(a) * Variables.annunr * drhonu) / a2;

      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgrho = dgrho + dgrho_phi;


      dgpres = (vpc.Variable_grhog(a) * deltag + vpc.Variable_grhor(a) * Variables.annur * deltar) / (a2 * 3.0) + vpc.Variable_grhonr(a) * Variables.annunr * dpnu / a2;

      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgpres = dgpres + dgprs_phi;

      if (Variables.initfl == 4) // Add a seed if desired.
         dgrho = dgrho + vpc.Variable_grhom(a) / a;



      dahdotdtau = -(dgrho + 3.0 * dgpres) * a;
      yprime[2] = dahdotdtau;
      // Force energy conservation.
      

      hdot = (2.0 * ak2 * eta + dgrho) / adotoa; // Eq(21a)

                                                           // 8*pi*G*(rho+P)*theta*a^2.
      dgtheta = vpc.Variable_grhom(a) * (Variables.omegac * thetac + Variables.omegab * thetab) / a + (4.0 / 3.0) * (vpc.Variable_grhog(a) * thetag + Variables.annur * vpc.Variable_grhor(a) * thetar) / a2 + Variables.annunr * vpc.Variable_grhonr(a) * ak * fnu / a2;

      // Variables.omegac,thetac,Variables.omegab,thetab,a,vpc.Variable_grhog(a), thetag,Variables.annur,vpc.Variable_grhor(a),thetar,a2,Variables.annunr,vpc.Variable_grhonr(a), ak,a2

      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgtheta = dgtheta + dgtheta_phi;

      etadot = 0.5 * dgtheta / ak2; //Eq(21b)
      if(fabs(Variables.omegak)>0.0001)
         etadot=0.5*(dgtheta + Variables.curv*hdot)/ak2;
      
     

      yprime[3] = etadot;


      if(fabs(Variables.omegak)<0.0001)
      {
         if (tcpa != 1) // First order approximation for photon shear
            shearg = (1.0 / (opac1 * 9.0)) * ((8.0 / 3.0) * thetag + (4.0 / 3.0) * hdot + 8.0 * etadot);
         }

      dgshear = (4.0 / 3.0) * (vpc.Variable_grhog(a) * shearg + Variables.annur * vpc.Variable_grhor(a) * shearr) / a2 + Variables.annunr * vpc.Variable_grhonr(a) * shearnu / a2; // 8*pi*G*(rho+P)*sigma*a**2.


      /*********************************************************************************************************************************/
      /*	  CDM equations of motion.																																		*/
      /*********************************************************************************************************************************/

      deltacdot = -thetac - 0.5 * hdot; // Ma & Bertschunger Eq. (29), (42) As not ineracting w=0 and delP=0
      yprime[4] = deltacdot;
      thetacdot = -adotoa * thetac; // Ma & Bertschunger Eq(29) As not ineracting w=0 and delP=0 K=0 as flat
      yprime[5] = thetacdot;



      /*********************************************************************************************************************************/
      /*  Photon equations of motion (total intensity).																											*/
      /*********************************************************************************************************************************/
      deltagdot = -(4.0 / 3.0) * (thetag + 0.5 * hdot); // Ma & Bertschunger Eq(63)
      yprime[8] = deltagdot;                            // As deltag in involved in baryon calculation so it is calculated here

      /*********************************************************************************************************************************/
      /*	  Baryon equations of motion.																																	*/
      /*********************************************************************************************************************************/
      deltabdot = -thetab - 0.5 * hdot; // Ma & Bertschunger Eq(66)
      yprime[6] = deltabdot;



      drag = opac1 * (thetag - thetab);

      if(fabs(Variables.omegak)>0.0001)
      {
         if( (tcp>0.999) && (fabs(Variables.omegak)>0.0001) )
            tcpb = 1;
         else 
            tcpb = 0;   
         }

      if (tcpb == 1)                                                                // If baryons and photons are uncoupled.
      {
         thetabdot = -adotoa * thetab + ak2 * cs2 * deltab + pb43 * drag; // Ma & Bertschunger Eq(66)
         //printf("Here IN TCPB");
         }
      else //If baryons and photons are tightly coupled. Zeroth-order approximation to baryon velocity.
      {
         thetabdot = (-adotoa * thetab + ak2 * cs2 * deltab + ak2 * pb43 * (0.25 * deltag - shearg)) / (1.0 + pb43); //
                                                                                                                                 //            printf("Line 505 : %e    %e    %e    %e    %e    %e    %e  %e\n",adotoa,thetab,Variables.ak2,cs2,deltab,pb43,deltag,shearg);
         if (Variables.dimflag == 1)
         {
            if (Variables.amnu == 0.0)
            {
               rhonudot = 0.0;
               shearnudot = 0.0;
            }
            else
            {
               Neutrino.nuder(a, rhonu, ysanpass, adotoa, y + Variables.iq2 - 1, yprime + Variables.iq2 - 1);
               rhonudot = ysanpass[0];
               shearnudot = ysanpass[1];
            }

            sgrhooa2 = sqrt(grho / a2);

            rgrho = sqrt(vpc.Variable_grhom(a) * vpc.Variable_Lambda(a) * vpc.Variable_Lambda(a) * Variables.omegav * vpc.Variable_grhom(a) * Variables.omegav + vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) / (a2 * a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / (a2 * a2));
            adotdota = sgrhooa2 / rgrho * (vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) * (-3 / (2 * a)) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) * (-2 / a2) + vpc.Variable_grhor(a) * Variables.annunr * rhonudot / (2 * a2 * adotoa)) + 2 * grho;

            adotdota = adotdota / 3.0;
         }

         else
         {
            adotdota = 0.5 * (adotoa * adotoa - gpres);
         }

         //  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
         slip = (2.0 * pb43 / (1.0 + pb43) + dlxedla) * adotoa * (thetab - thetag) + (1.0 / opac1) * (-adotdota * thetab - adotoa * ak2 * 0.5 * deltag + ak2 * (cs2 * deltabdot - 0.25 * deltagdot)) / (1.0 + pb43);
         //slip = (2.0 * pb43 / (1.0 + pb43)) * adotoa * (thetab - thetag) + (1.0 / opac1) * (-adotdota * thetab - adotoa * ak2 * 0.5 * deltag + ak2 * (cs2 * deltabdot - 0.25 * deltagdot)) / (1.0 + pb43);         
         //  First-order approximation to baryon velocity.

         thetabdot = thetabdot + pb43 / (1.0 + pb43) * slip;
      }


      yprime[7] = thetabdot;


      /*********************************************************************************************************************************/
      /*	  Photon equations of motion (polarization).  																											*/
      /*********************************************************************************************************************************/
      thetagdot = (-thetabdot - adotoa * thetab + ak2 * cs2 * deltab) / pb43 + ak2 * (0.25 * deltag -  Curvature.b(2,ak2) * shearg);
      yprime[9] = thetagdot; // Ma & Bertschunger Eq(75)

      if(Variables.VPCflag == 1)
      {
         yprime[4] = yprime[4] - vpc.Variable_G_c(a) * y[4];
         yprime[5] = yprime[5] - vpc.Variable_G_c(a) * y[5];
         yprime[6] = yprime[6] - vpc.Variable_G_c(a) * y[6];
         yprime[7] = yprime[7] - vpc.Variable_G_c(a) * y[7];
         yprime[8] = yprime[8] - vpc.Variable_G_c(a) * y[8];
         yprime[9] = yprime[9] - vpc.Variable_G_c(a) * y[9];
         }


      if(fabs(Variables.omegak)>0.0001)
      {
         if(tcp>0.999)
            tcpa = 1;
         else
            tcpa = 0;
         }

      //printf("\n%e %e %e %e %e %e %e %e",a,eta,deltac,thetac,deltab,thetab,deltag,thetag);
      //printf("\n%e %e",dgrho,dgpres);
      //printf("\n%e %e %e %e %e",yprime[1],yprime[2],yprime[3],yprime[4],yprime[5]);
      //printf("\n%e %e %e %e %e\n",yprime[6],yprime[7],yprime[8],yprime[9],yprime[10]);      
      //exit(1);

      if (tcpa == 1) // Treat baryons and photons as uncoupled.
      {
         yprime[10] = (8.0 / 15.0) * Curvature.b(2,ak2) *  thetag - 0.6 * ak * Curvature.b(3,ak2) *  y[11] - opac1 * y[10] + (4.0 / 15.0) * Curvature.b(2,ak2) * hdot + (8.0 / 5.0) * Curvature.b(2,ak2) * etadot * Variables.aux + 0.1 * opac1 * polter;
         // Polarization equations for l = 0, 1, 2.

         yprime[9 + Variables.lmaxg] = -ak * y[10 + Variables.lmaxg] - opac1 * y[9 + Variables.lmaxg] + 0.5 * opac1 * polter;
         yprime[10 + Variables.lmaxg] = (ak / 3.0) * (y[9 + Variables.lmaxg] - 2.0 * Curvature.b(2,ak2) * y[11 + Variables.lmaxg]) - opac1 * y[10 + Variables.lmaxg];
         yprime[11 + Variables.lmaxg] = ak * (0.4 * Curvature.b(2,ak2) * y[10 + Variables.lmaxg] - 0.6 * Curvature.b(3,ak2) * y[12 + Variables.lmaxg]) - opac1 * y[11 + Variables.lmaxg] + 0.1 * opac1 * polter;

         if(Variables.VPCflag == 1)
         {
            yprime[10] = yprime[10] - vpc.Variable_G_c(a) * y[10];
            yprime[9 + Variables.lmaxg] = yprime[9 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg];
            yprime[10 + Variables.lmaxg] = yprime[10 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + Variables.lmaxg];
            yprime[11 + Variables.lmaxg] = yprime[11 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + Variables.lmaxg];
            }

         for (int l = 3; l <= Variables.lmaxg - 1; l++)
         {
            yprime[8 + l] = ak * (1.0 / (2 * l + 1.0)) * (l * Curvature.b(l,ak2) * y[7 + l] - (l + 1) * Curvature.b(l+1,ak2) * y[9 + l]) - opac1 * y[8 + l];
            yprime[9 + Variables.lmaxg + l] = ak * (1.0 / (2 * l + 1.0)) * (l * Curvature.b(l,ak2) * y[8 + Variables.lmaxg + l] - (l + 1) * Curvature.b(l+1,ak2) * y[10 + Variables.lmaxg + l]) - opac1 * y[9 + Variables.lmaxg + l];
            if(Variables.VPCflag == 1)
            {
               yprime[8 + l] = yprime[8 + l] - vpc.Variable_G_c(a) * y[8 + l];
               yprime[9 + Variables.lmaxg + l] = yprime[9 + Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg + l];
               }
         }
      }
      else // Treat baryons and photons as tightly coupled (with no polarization).
      {
         yprime[10] = 0.0;
         yprime[9 + Variables.lmaxg] = 0.0;
         yprime[10 + Variables.lmaxg] = 0.0;
         yprime[11 + Variables.lmaxg] = 0.0;

         if(Variables.VPCflag == 1)
         {
            yprime[10] = yprime[10] - vpc.Variable_G_c(a) * y[10];
            yprime[9 + Variables.lmaxg] = yprime[9 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg];
            yprime[10 + Variables.lmaxg] = yprime[10 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + Variables.lmaxg];
            yprime[11 + Variables.lmaxg] = yprime[11 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + Variables.lmaxg];
            }

         for (int l = 3; l <= Variables.lmaxg - 1; l++)
         {
            yprime[8 + l] = 0.0;
            yprime[9 + Variables.lmaxg + l] = 0.0;
            if(Variables.VPCflag == 1)
            {
               yprime[8 + l] = yprime[8 + l] - vpc.Variable_G_c(a) * y[8 + l];
               yprime[9 + Variables.lmaxg + l] = yprime[9 + Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg + l];
               }            
         }
      } // Truncate moment expansion.

      //printf("\n--------- B:%e %e %e %e %e %e",dgshear,a,vpc.Variable_grhog(a),shearg,shearr,yprime[10]);

      //printf("\nYP %e %e %e %e",tau,yprime[7],yprime[8],yprime[9]);//,yprime[10],yprime[11]);
      yprime[8 + Variables.lmaxg] = ak * Curvature.b(Variables.lmaxg, ak2) * y[7 + Variables.lmaxg] - (Variables.lmaxg + 1) / tau * y[8 + Variables.lmaxg] - opac1 * y[8 + Variables.lmaxg];
      yprime[9 + 2 * Variables.lmaxg] = ak * Curvature.b(Variables.lmaxg, ak2) * y[8 + 2 * Variables.lmaxg] - (Variables.lmaxg + 1) / tau * y[9 + 2 * Variables.lmaxg] - opac1 * y[9 + 2 * Variables.lmaxg];

      /*********************************************************************************************************************************/
      /* 	 Massless Neutrino equations of motion.																												*/
      /*********************************************************************************************************************************/
      deltardot = -(4.0 / 3.0) * (thetar + 0.5 * hdot); // Ma & Bertschunger Eq(49a)
      yprime[10 + 2 * Variables.lmaxg] = deltardot;

      thetardot = ak2 * (0.25 * deltar - Curvature.b(2, ak2) *shearr); // Ma & Bertschunger Eq(49b)
      yprime[11 + 2 * Variables.lmaxg] = thetardot;

      Fnu2dot = 8.0 / 15.0 * Curvature.b(2,ak2) * thetar - 0.6 * ak * y[13 + 2 * Variables.lmaxg] * Curvature.b(3,ak2) + 4.0 / 15.0 * hdot * Curvature.b(2,ak2) + 8.0 / 5.0 * etadot * Variables.aux * Curvature.b(2, ak2);
      yprime[12 + 2 * Variables.lmaxg] = Fnu2dot; // Ma & Bertschunger Eq(49c)

      if(Variables.VPCflag == 1)
      {
         yprime[8 + Variables.lmaxg] = yprime[8 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[8 + Variables.lmaxg];
         yprime[9 + 2 * Variables.lmaxg] = yprime[9 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + 2 * Variables.lmaxg];
         yprime[10 + 2 * Variables.lmaxg] = yprime[10 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg];
         yprime[11 + 2 * Variables.lmaxg] = yprime[11 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + 2 * Variables.lmaxg];
         yprime[12 + 2 * Variables.lmaxg] = yprime[12 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[12 + 2 * Variables.lmaxg];         
         }


      for (int l = 3; l <= Variables.lmaxnr - 1; l++)
      {
         yprime[10 + 2 * Variables.lmaxg + l] = ak * (l * Curvature.b(l,ak2) * y[9 + 2 * Variables.lmaxg + l] - (l + 1) * Curvature.b(l + 1, ak2) * y[11 + 2 * Variables.lmaxg + l]) / (2 * l + 1); // Ma & Bertschunger Eq(49d)

         if(Variables.VPCflag == 1)
            yprime[10 + 2 * Variables.lmaxg + l] = yprime[10 + 2 * Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg + l];
         }

      // Truncate moment expansion. Eq(51)
      yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] = ak * Curvature.b(Variables.lmaxg, ak2) * y[9 + 2 * Variables.lmaxg + Variables.lmaxnr] - (Variables.lmaxnr + 1) / tau * y[10 + 2 * Variables.lmaxg + Variables.lmaxnr];
      
      if(Variables.VPCflag == 1)
         yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] = yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg + Variables.lmaxnr];

      /*********************************************************************************************************************************/
      /*	  Equations for dark energy perturbations.																												*/
      /*********************************************************************************************************************************/
      if ((Variables.ndyn == 1) || (Variables.ndyn == 2))
      {
         Darkenergy.dyn_phi(a, hdot, grho, gpres, rc_phi, rc_psi, ysanpass); // Dark Energy perturbation
         yprime[n - 1] = ysanpass[0];                                        //rc_dphi;
         yprime[n] = ysanpass[1];                                            //rc_dpsi;
         if(Variables.VPCflag == 1)
         {
            yprime[n-1] = yprime[n-1] - vpc.Variable_G_c(a) * y[n-1];
            yprime[n] = yprime[n] - vpc.Variable_G_c(a) * y[n];            
            }
      }

      //for(int ig=1;ig<40;ig++)
      //   printf("\n -- %d %e %e",ig,y[ig],yprime[ig]);   

      //if(Variables.VPCflag == 1)
      //{
       //  for(int i=0;i<n;i++)
       //     yprime[i] = yprime[i] * vpc.Variable_C(a);
       //  }

      //printf("\n%e %e %e %e %e %e",a,tcp,tcp1,tcp2,y[9],y[8]);

      if (Variables.nqmax == 0)
         return; //Return if no massive Neutrinos are present
                 //         printf("Line 610\n");
      /*********************************************************************************************************************************/
      /*  Massive Neutrino equations of motion.																														*/
      /*********************************************************************************************************************************/

      dq = 1.0;

      for (int i = 1; i <= Variables.nqmax; i++)
      {
         q = i * dq - 0.5;
         aq = a * Variables.amnu / q;
         v = 1.0 / sqrt(1.0 + aq * aq);
         akv[i] = ak * v;
      }

      for (int i = 1; i <= Variables.nqmax; i++)
      {
         ind = Variables.iq0 + i - 1;
         yprime[ind] = -akv[i] * y[ind + Variables.nqmax] + hdot * GlobalArray.dlfdlq[i] / 6.0; // Ma & Bertschunger Eq(56a)
         if(Variables.VPCflag == 1)  
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind];

         ind = Variables.iq1 + i - 1;
         yprime[ind] = akv[i] * (y[ind - Variables.nqmax] - 2 * Curvature.b(2,ak2) * y[ind + Variables.nqmax]) / 3; // Ma & Bertschunger Eq(56b)
         if(Variables.VPCflag == 1)   
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 

         ind = Variables.iq2 + i - 1;                                                          // Ma & Bertschunger Eq(56c)
         yprime[ind] = akv[i] * (2 * Curvature.b(2,ak2) * y[ind - Variables.nqmax] - 3 * Curvature.b(3,ak2) * y[ind + Variables.nqmax]) / 5 - (hdot / 15.0 + 2.0 / 5.0 * etadot * Variables.aux ) * Curvature.b(2,ak2) * GlobalArray.dlfdlq[i];
         if(Variables.VPCflag == 1)   
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind];
         
         ind = 10 + 2 * Variables.lmaxg + Variables.lmaxnr + i + Variables.lmaxnu * Variables.nqmax;          // Truncate moment expansion.
         yprime[ind] = akv[i] * y[ind - Variables.nqmax]*Curvature.b(i,ak2) - (Variables.lmaxnu + 1) / tau * y[ind]; // Ma & Bertschunger Eq(57c) Eq(58) Combine
         if(Variables.VPCflag == 1)   
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind];
      }

      for (int l = 3; l <= Variables.lmaxnu - 1; l++)
      {
         for (int i = 1; i <= Variables.nqmax; i++)
         {
            ind = 10 + 2 * Variables.lmaxg + Variables.lmaxnr + i + l * Variables.nqmax;
            yprime[ind] = akv[i] * (1.0 / (2 * l + 1.0)) * (l * Curvature.b(l,ak2) *y[ind - Variables.nqmax] - (l + 1) * Curvature.b(l + 1,ak2) * y[ind + Variables.nqmax]); // Ma & Bertschunger Eq(57c)
            if(Variables.VPCflag == 1)
            {   
               yprime[ind] = (yprime[ind] - vpc.Variable_G_c(a) * y[ind]);
            //   yprime[ind] = vpc.Variable_C(a) * yprime[ind];
               }
         }
      }

   }



   /************************************************************************************************************************************/
   /* 	This will calculate the derivatives of the scaler perturbation terms. This part is totally based on Ma & Bertschunger			*/
   /************************************************************************************************************************************/

   void fderivs(int n, double x, double y[], double yprime[])
   { // n store total number of values in y
      double ep0 = .01;
      double xsanpass[3], ep, tau, epstc, dpnu, grho, gpres;
      double deltac, thetac, deltab, thetab, deltag, thetag, shearg;
      double polter, deltar, thetar, shearr, a2, cs2, dlxedla; //,opec
      double a, eta, rc_phi, rc_psi, photbar, pb43, ysanpass[2];
      double tcp, tcp1, tcp2, tcp3, drag, deltabdot;
      double rgrho, slip, pnu, fnu, shearnu;
      double dgrho, dgpres, dgtheta, etadot;
      double deltacdot, thetacdot, deltagdot, thetabdot;
      double rhonudot, shearnudot, adotdota, weos;
      double thetagdot, deltardot, thetardot, Fnu2dot;
      double dgprs_phi, dahdotdtau, sgrhooa2;
      double opac1, dgrho_phi, dgtheta_phi, dq, q, v, aq;
      double akv[nqmax0 + 1]; //,rc_dphi,rc_dpsi;
      double zsanpass[4];

      int ind, tcpa, tcpb;
      VPC vpc;


      neutrino Neutrino;
      darkenergy Darkenergy;
      numericx Numericx;
      recombination Recombination;
      if (Variables.ak > Variables.epsw) // ep is used to stop the tight coupling approximation.
         ep = ep0;
      else
         ep = 0.5 * ep0;

      tau = x; // time of the derivative
      
      //if(testVariable%100 == 0)
      //   printf("\n\t\t %d %e", testVariable++,tau);

      a = y[1];
      eta = y[3];

      deltac = y[4]; //  CDM.
      thetac = y[5];

      deltab = y[6]; // Baryons.
      thetab = y[7];

      deltag = y[8]; //  Photons.
      thetag = y[9];
      shearg = y[10] / 2.0;

      polter = y[10] + y[9 + Variables.lmaxg] + y[11 + Variables.lmaxg]; //  Polarization term. Ma Eq(60). F2+G0+G2

      deltar = y[10 + 2 * Variables.lmaxg]; //  Massless neutrinos.
      thetar = y[11 + 2 * Variables.lmaxg];
      shearr = y[12 + 2 * Variables.lmaxg] / 2.0;
      
      //printf("\n\nThis tau= %e %e\n",tau, dlxedla);

      a2 = a * a;
      Recombination.thermo(tau, xsanpass); // Get the Value of dlxedla

      cs2 = xsanpass[0];
      opac1 = xsanpass[1];
      dlxedla = xsanpass[2];

      //printf("This tau= %e %e\n\n",tau, dlxedla);

      photbar = vpc.Variable_grhog(a) / (vpc.Variable_grhom(a) * Variables.omegab * a); // Photon mass density over baryon mass density.
      pb43 = (4.0 / 3.0) * photbar;

      tcp = 0.0; // Tight Coupling parameters
      tcp1 = Variables.ak / opac1;
      tcp2 = 1.0 / (opac1 * tau);
      tcp3 = Variables.ak * tau;

      tcpa = 0;
      tcpb = 0;

      epstc = Numericx.max(tcp1, tcp2);

      if ((epstc > 5.0e-3) && (a > 1.35e-5))
         tcpa = 1;

      epstc = epstc / (1.0 + pb43);

      if ((epstc > 5.0e-3) && (a > 1.35e-4))
         tcpb = 1;

      /*********************************************************************************************************************************/
      /*  Compute expansion rate (a dot = sqrt(8 pi G rho / 3)).																								*/
      /*********************************************************************************************************************************/

      if (Variables.amnu == 0.0) // No massive Neutrino, so no perturbation
      {
         rhonu = 1.0;
         pnu = 1.0 / 3.0;
         drhonu = 0.0;
         fnu = 0.0;
         dpnu = 0.0;
         shearnu = 0.0;
      }
      else // If massive Neutrino present then perturbe
      {
         Neutrino.nu1(a, ysanpass);
         rhonu = ysanpass[0];
         pnu = ysanpass[1];
         Neutrino.nu2(a, zsanpass, y + Variables.iq0, y + Variables.iq1, y + Variables.iq2);
         drhonu = zsanpass[0];
         fnu = zsanpass[1];
         dpnu = zsanpass[2];
         shearnu = zsanpass[3];
      }

      Variables.omegavdyn = Variables.omegav * Darkenergy.dynrho(a);
      
      if(Variables.VPCflag == 1)  
         Variables.omegavdyn = Variables.omegavdyn * vpc.Variable_Lambda(a);


      weos = Darkenergy.wdyn_func(a);

      //      printf("Line 414\n");																								//  8*pi*G*rho*a^2 and 8*pi*G*P*a^2.

      grho = vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + Variables.annunr * rhonu * vpc.Variable_grhonr(a)) / a2 + vpc.Variable_grhom(a) * Variables.omegavdyn * a2;

      if (Variables.dimflag == 1)
         grho = pow((sqrt(grho) + sqrt(Variables.omegav * vpc.Variable_grhom(a) * a2)), 2); // energy density for  5 dimensions

      adotoa = sqrt(grho / 3.0);
      yprime[1] = adotoa * a;

      gpres = ((vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur) / 3.0 + vpc.Variable_grhonr(a) * Variables.annunr * pnu) / a2; // 8*pi*G*rho*a^2 and 8*pi*G*P*a^2.
      gpres = gpres + weos * vpc.Variable_grhom(a) * Variables.omegavdyn * a2;

      if (Variables.ndyn == 1 || Variables.ndyn == 2) // dark energy perturbations
      {
         rc_phi = y[n - 1];
         rc_psi = y[n];

         Darkenergy.dyn_nrg(a, grho, gpres, rc_phi, rc_psi, xsanpass);
         dgrho_phi = xsanpass[0];
         dgprs_phi = xsanpass[1];
         dgtheta_phi = xsanpass[2];
      }
      // Evaluate metric and massive neutrino perturbations.
      //  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
      dgrho = vpc.Variable_grhom(a) * (Variables.omegac * deltac + Variables.omegab * deltab) / a + (vpc.Variable_grhog(a) * deltag + vpc.Variable_grhor(a) * Variables.annur * deltar + vpc.Variable_grhonr(a) * Variables.annunr * drhonu) / a2;

      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgrho = dgrho + dgrho_phi;

      dgpres = (vpc.Variable_grhog(a) * deltag + vpc.Variable_grhor(a) * Variables.annur * deltar) / (a2 * 3.0) + vpc.Variable_grhonr(a) * Variables.annunr * dpnu / a2;

      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgpres = dgpres + dgprs_phi;

      if (Variables.initfl == 4) // Add a seed if desired.
         dgrho = dgrho + vpc.Variable_grhom(a) / a;

      dahdotdtau = -(dgrho + 3.0 * dgpres) * a;
      yprime[2] = dahdotdtau;
      // Force energy conservation.
      hdot = (2.0 * Variables.ak2 * eta + dgrho) / adotoa; // Eq(21a)
                                                           // 8*pi*G*(rho+P)*theta*a^2.
      dgtheta = vpc.Variable_grhom(a) * (Variables.omegac * thetac + Variables.omegab * thetab) / a + (4.0 / 3.0) * (vpc.Variable_grhog(a) * thetag + Variables.annur * vpc.Variable_grhor(a) * thetar) / a2 + Variables.annunr * vpc.Variable_grhonr(a) * Variables.ak * fnu / a2;


      if (Variables.ndyn == 1 || Variables.ndyn == 2)
         dgtheta = dgtheta + dgtheta_phi;

      etadot = 0.5 * dgtheta / Variables.ak2; //Eq(21b)
      yprime[3] = etadot;

      if (tcpa != 1) // First order approximation for photon shear
         shearg = (1.0 / (opac1 * 9.0)) * ((8.0 / 3.0) * thetag + (4.0 / 3.0) * hdot + 8.0 * etadot);

      dgshear = (4.0 / 3.0) * (vpc.Variable_grhog(a) * shearg + Variables.annur * vpc.Variable_grhor(a) * shearr) / a2 + Variables.annunr * vpc.Variable_grhonr(a) * shearnu / a2; // 8*pi*G*(rho+P)*sigma*a**2.

      /*********************************************************************************************************************************/
      /*	  CDM equations of motion.																																		*/
      /*********************************************************************************************************************************/

      deltacdot = -thetac - 0.5 * hdot; // Ma & Bertschunger Eq. (29), (42) As not ineracting w=0 and delP=0
      yprime[4] = deltacdot;
      thetacdot = -adotoa * thetac; // Ma & Bertschunger Eq(29) As not ineracting w=0 and delP=0 K=0 as flat
      yprime[5] = thetacdot;

      /*********************************************************************************************************************************/
      /*  Photon equations of motion (total intensity).																											*/
      /*********************************************************************************************************************************/
      deltagdot = -(4.0 / 3.0) * (thetag + 0.5 * hdot); // Ma & Bertschunger Eq(63)
      yprime[8] = deltagdot;                            // As deltag in involved in baryon calculation so it is calculated here

      /*********************************************************************************************************************************/
      /*	  Baryon equations of motion.																																	*/
      /*********************************************************************************************************************************/
      deltabdot = -thetab - 0.5 * hdot; // Ma & Bertschunger Eq(66)
      yprime[6] = deltabdot;

      drag = opac1 * (thetag - thetab);


      if (tcpb == 1)                                                                // If baryons and photons are uncoupled.
         thetabdot = -adotoa * thetab + Variables.ak2 * cs2 * deltab + pb43 * drag; // Ma & Bertschunger Eq(66)

      else //If baryons and photons are tightly coupled. Zeroth-order approximation to baryon velocity.
      {
         thetabdot = (-adotoa * thetab + Variables.ak2 * cs2 * deltab + Variables.ak2 * pb43 * (0.25 * deltag - shearg)) / (1.0 + pb43); //
                                                                                                                                         //            printf("Line 505 : %e    %e    %e    %e    %e    %e    %e  %e\n",adotoa,thetab,Variables.ak2,cs2,deltab,pb43,deltag,shearg);
         if (Variables.dimflag == 1)
         {
            if (Variables.amnu == 0.0)
            {
               rhonudot = 0.0;
               shearnudot = 0.0;
            }
            else
            {
               Neutrino.nuder(a, rhonu, ysanpass, adotoa, y + Variables.iq2 - 1, yprime + Variables.iq2 - 1);
               rhonudot = ysanpass[0];
               shearnudot = ysanpass[1];
            }

            sgrhooa2 = sqrt(grho / a2);

            rgrho = sqrt(vpc.Variable_grhom(a) * Variables.omegav * vpc.Variable_grhom(a) * Variables.omegav  * vpc.Variable_Lambda(a)  * vpc.Variable_Lambda(a) + vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) / (a2 * a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / (a2 * a2));
            adotdota = sgrhooa2 / rgrho * (vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) * (-3 / (2 * a)) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) * (-2 / a2) + vpc.Variable_grhor(a) * Variables.annunr * rhonudot / (2 * a2 * adotoa)) + 2 * grho;

            adotdota = adotdota / 3.0;
         }

         else
         {
            adotdota = 0.5 * (adotoa * adotoa - gpres);
         }


         //  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
         slip = (2.0 * pb43 / (1.0 + pb43) + dlxedla) * adotoa * (thetab - thetag) + (1.0 / opac1) * (-adotdota * thetab - adotoa * Variables.ak2 * 0.5 * deltag + Variables.ak2 * (cs2 * deltabdot - 0.25 * deltagdot)) / (1.0 + pb43);
         //  First-order approximation to baryon velocity.

         thetabdot = thetabdot + pb43 / (1.0 + pb43) * slip;
      }

      yprime[7] = thetabdot;

      /*********************************************************************************************************************************/
      /*	  Photon equations of motion (polarization).  																											*/
      /*********************************************************************************************************************************/
      thetagdot = (-thetabdot - adotoa * thetab + Variables.ak2 * cs2 * deltab) / pb43 + Variables.ak2 * (0.25 * deltag - shearg);
      yprime[9] = thetagdot; // Ma & Bertschunger Eq(75)

      if(Variables.VPCflag == 1)
      {
         yprime[4] = yprime[4] - vpc.Variable_G_c(a) * y[4];
         yprime[5] = yprime[5] - vpc.Variable_G_c(a) * y[5];
         yprime[6] = yprime[6] - vpc.Variable_G_c(a) * y[6];
         yprime[7] = yprime[7] - vpc.Variable_G_c(a) * y[7];
         yprime[8] = yprime[8] - vpc.Variable_G_c(a) * y[8];
         yprime[9] = yprime[9] - vpc.Variable_G_c(a) * y[9];
         }

      if (tcpa == 1) // Treat baryons and photons as uncoupled.
      {
         yprime[10] = (8.0 / 15.0) * thetag - 0.6 * Variables.ak * y[11] - opac1 * y[10] + (4.0 / 15.0) * hdot + (8.0 / 5.0) * etadot + 0.1 * opac1 * polter;
         // Polarization equations for l = 0, 1, 2.
         yprime[9 + Variables.lmaxg] = -Variables.ak * y[10 + Variables.lmaxg] - opac1 * y[9 + Variables.lmaxg] + 0.5 * opac1 * polter;
         yprime[10 + Variables.lmaxg] = (Variables.ak / 3.0) * (y[9 + Variables.lmaxg] - 2.0 * y[11 + Variables.lmaxg]) - opac1 * y[10 + Variables.lmaxg];
         yprime[11 + Variables.lmaxg] = Variables.ak * (0.4 * y[10 + Variables.lmaxg] - 0.6 * y[12 + Variables.lmaxg]) - opac1 * y[11 + Variables.lmaxg] + 0.1 * opac1 * polter;

         if(Variables.VPCflag == 1)
         {
            yprime[10] = yprime[10] - vpc.Variable_G_c(a) * y[10];
            yprime[9 + Variables.lmaxg] = yprime[9 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg];
            yprime[10 + Variables.lmaxg] = yprime[10 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + Variables.lmaxg];
            yprime[11 + Variables.lmaxg] = yprime[11 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + Variables.lmaxg];
            }

         for (int l = 3; l <= Variables.lmaxg - 1; l++)
         {
            yprime[8 + l] = Variables.ak * (1.0 / (2 * l + 1.0)) * (l * y[7 + l] - (l + 1) * y[9 + l]) - opac1 * y[8 + l];
            yprime[9 + Variables.lmaxg + l] = Variables.ak * (1.0 / (2 * l + 1.0)) * (l * y[8 + Variables.lmaxg + l] - (l + 1) * y[10 + Variables.lmaxg + l]) - opac1 * y[9 + Variables.lmaxg + l];
            if(Variables.VPCflag == 1)
            {
               yprime[8 + l] = yprime[8 + l] - vpc.Variable_G_c(a) * y[8 + l];
               yprime[9 + Variables.lmaxg + l] = yprime[9 + Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg + l];
               }
         }
      }
      else // Treat baryons and photons as tightly coupled (with no polarization).
      {
         yprime[10] = 0.0;
         yprime[9 + Variables.lmaxg] = 0.0;
         yprime[10 + Variables.lmaxg] = 0.0;
         yprime[11 + Variables.lmaxg] = 0.0;

         if(Variables.VPCflag == 1)
            {
               yprime[10] = yprime[10] - vpc.Variable_G_c(a) * y[10];
               yprime[9 + Variables.lmaxg] = yprime[9 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg];
               yprime[10 + Variables.lmaxg] = yprime[10 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + Variables.lmaxg];
               yprime[11 + Variables.lmaxg] = yprime[11 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + Variables.lmaxg];
               }

         for (int l = 3; l <= Variables.lmaxg - 1; l++)
         {
            yprime[8 + l] = 0.0;
            yprime[9 + Variables.lmaxg + l] = 0.0;
            
            if(Variables.VPCflag == 1)
            {
               yprime[8 + l] = yprime[8 + l] - vpc.Variable_G_c(a) * y[8 + l];
               yprime[9 + Variables.lmaxg + l] = yprime[9 + Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[9 + Variables.lmaxg + l];
               }              
         }
      } // Truncate moment expansion.

      yprime[8 + Variables.lmaxg] = Variables.ak * y[7 + Variables.lmaxg] - (Variables.lmaxg + 1) / tau * y[8 + Variables.lmaxg] - opac1 * y[8 + Variables.lmaxg];
      yprime[9 + 2 * Variables.lmaxg] = Variables.ak * y[8 + 2 * Variables.lmaxg] - (Variables.lmaxg + 1) / tau * y[9 + 2 * Variables.lmaxg] - opac1 * y[9 + 2 * Variables.lmaxg];

      /*********************************************************************************************************************************/
      /* 	 Massless Neutrino equations of motion.																												*/
      /*********************************************************************************************************************************/
      deltardot = -(4.0 / 3.0) * (thetar + 0.5 * hdot); // Ma & Bertschunger Eq(49a)
      yprime[10 + 2 * Variables.lmaxg] = deltardot;

      thetardot = Variables.ak2 * (0.25 * deltar - shearr); // Ma & Bertschunger Eq(49b)
      yprime[11 + 2 * Variables.lmaxg] = thetardot;

      Fnu2dot = 8.0 / 15.0 * thetar - 0.6 * Variables.ak * y[13 + 2 * Variables.lmaxg] + 4.0 / 15.0 * hdot + 8.0 / 5.0 * etadot;
      yprime[12 + 2 * Variables.lmaxg] = Fnu2dot; // Ma & Bertschunger Eq(49c)

      if(Variables.VPCflag == 1)
      {
         yprime[8 + Variables.lmaxg] = yprime[8 + Variables.lmaxg] - vpc.Variable_G_c(a) * y[8 + Variables.lmaxg];
         yprime[9 + 2 * Variables.lmaxg] = yprime[9 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[9 + 2 * Variables.lmaxg];
         yprime[10 + 2 * Variables.lmaxg] = yprime[10 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg];
         yprime[11 + 2 * Variables.lmaxg] = yprime[11 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[11 + 2 * Variables.lmaxg];
         yprime[12 + 2 * Variables.lmaxg] = yprime[12 + 2 * Variables.lmaxg] - vpc.Variable_G_c(a) * y[12 + 2 * Variables.lmaxg];         
         }  


      for (int l = 3; l <= Variables.lmaxnr - 1; l++)
      {
         yprime[10 + 2 * Variables.lmaxg + l] = Variables.ak * (l * y[9 + 2 * Variables.lmaxg + l] - (l + 1) * y[11 + 2 * Variables.lmaxg + l]) / (2 * l + 1); // Ma & Bertschunger Eq(49d)
         
         if(Variables.VPCflag == 1)
            yprime[10 + 2 * Variables.lmaxg + l] = yprime[10 + 2 * Variables.lmaxg + l] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg + l];         
         }
      // Truncate moment expansion. Eq(51)
      yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] = Variables.ak * y[9 + 2 * Variables.lmaxg + Variables.lmaxnr] - (Variables.lmaxnr + 1) / tau * y[10 + 2 * Variables.lmaxg + Variables.lmaxnr];
      if(Variables.VPCflag == 1)
         yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] = yprime[10 + 2 * Variables.lmaxg + Variables.lmaxnr] - vpc.Variable_G_c(a) * y[10 + 2 * Variables.lmaxg + Variables.lmaxnr];


      /*********************************************************************************************************************************/
      /*	  Equations for dark energy perturbations.																												*/
      /*********************************************************************************************************************************/
      if ((Variables.ndyn == 1) || (Variables.ndyn == 2))
      {
         Darkenergy.dyn_phi(a, hdot, grho, gpres, rc_phi, rc_psi, ysanpass); // Dark Energy perturbation
         yprime[n - 1] = ysanpass[0];                                        //rc_dphi;
         yprime[n] = ysanpass[1];       
         
         if(Variables.VPCflag == 1)
         {
            yprime[n-1] = yprime[n-1] - vpc.Variable_G_c(a) * y[n-1];
            yprime[n] = yprime[n] - vpc.Variable_G_c(a) * y[n];            
            }                                              //rc_dpsi;
      }

      //for(int iii=1;iii<50;iii++)
      //   printf("\n %d %e %e",iii, y[iii],yprime[iii]);

      if (Variables.nqmax == 0)
         return; //Return if no massive Neutrinos are present
                 //         printf("Line 610\n");
      /*********************************************************************************************************************************/
      /*  Massive Neutrino equations of motion.																														*/
      /*********************************************************************************************************************************/

      dq = 1.0;

      for (int i = 1; i <= Variables.nqmax; i++)
      {
         q = i * dq - 0.5;
         aq = a * Variables.amnu / q;
         v = 1.0 / sqrt(1.0 + aq * aq);
         akv[i] = Variables.ak * v;
      }

      for (int i = 1; i <= Variables.nqmax; i++)
      {
         ind = Variables.iq0 + i - 1;
         yprime[ind] = -akv[i] * y[ind + Variables.nqmax] + hdot * GlobalArray.dlfdlq[i] / 6.0; // Ma & Bertschunger Eq(56a)
         if(Variables.VPCflag == 1)
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 

         ind = Variables.iq1 + i - 1;
         yprime[ind] = akv[i] * (y[ind - Variables.nqmax] - 2 * y[ind + Variables.nqmax]) / 3; // Ma & Bertschunger Eq(56b)
         if(Variables.VPCflag == 1)
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 

         ind = Variables.iq2 + i - 1;                                                          // Ma & Bertschunger Eq(56c)
         yprime[ind] = akv[i] * (2 * y[ind - Variables.nqmax] - 3 * y[ind + Variables.nqmax]) / 5 - (hdot / 15.0 + 2.0 / 5.0 * etadot) * GlobalArray.dlfdlq[i];
         if(Variables.VPCflag == 1)
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 

         ind = 10 + 2 * Variables.lmaxg + Variables.lmaxnr + i + Variables.lmaxnu * Variables.nqmax;
         // Truncate moment expansion.
         yprime[ind] = akv[i] * y[ind - Variables.nqmax] - (Variables.lmaxnu + 1) / tau * y[ind]; // Ma & Bertschunger Eq(57c) Eq(58) Combine
         if(Variables.VPCflag == 1)
            yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 
      }

      for (int l = 3; l <= Variables.lmaxnu - 1; l++)
      {
         for (int i = 1; i <= Variables.nqmax; i++)
         {
            ind = 10 + 2 * Variables.lmaxg + Variables.lmaxnr + i + l * Variables.nqmax;
            yprime[ind] = akv[i] * (1.0 / (2 * l + 1.0)) * (l * y[ind - Variables.nqmax] - (l + 1) * y[ind + Variables.nqmax]); // Ma & Bertschunger Eq(57c)
            if(Variables.VPCflag == 1)
               yprime[ind] = yprime[ind] - vpc.Variable_G_c(a) * y[ind]; 
         }
      }

   }





   void foutputCurved(int n, double y[], double yprime[], int j, double tau0, double tau, double zsanpass[])
   {
      double x, a, a2, eta, etadot, alpha, alphadot, deltab, thetab;
      double thetabdot, deltag, thetagdot, sheargdot, polter, coupldot;
      double polterddot, shearrdot, rhonudot, shearnudot, omegavdyn;
      double weos, grhodot, grho, sgrhooa2, rgrho, adotoadot;
      double dgsheardot, alphaddot, s1, s2, chi, chir, glens;
      double polterdot, ysanpass[2];
      double d1, dp1, dk1, phi1;

      double ak,ak2;
      double rsinh2;
      curvature Curvature = curvature(-Variables.omegak);
      VPC vpc;

      ak = Variables.ak;
      ak2 = Variables.ak2;

      if(fabs(Variables.omegak)>0.0001)
      {
         ak = sqrt(ak2 - Variables.curv);
         ak2 = ak*ak;
      }

      chi=(tau0-tau)/Variables.r;
      rsinh2=(Variables.r*Curvature.sinhK(chi));
      rsinh2 = rsinh2 * rsinh2;     


      //printf("\n\nRzinh: %e",rsinh2);
      //exit(1);   

      darkenergy Darkenergy;
      neutrino Neutrino;
      //curvature Curvature(-Variables.omegak);
      d1 = 0.0;
      dp1 = 0.0;
      dk1 = 0.0;

      //printf("\n\n%e ",dgshear);
      //exit(1);

      x = ak * (tau0 - tau);
      a = y[1];
      a2 = a * a;

      eta = y[3];
      etadot = yprime[3];

      alpha = (hdot + 6 * etadot*Variables.aux) / (2.0 * ak2);                          // M.Zaldarriaga's Ph.D. Thesis P.63
      alphadot = -3 * dgshear / (2.0 * ak2 *Curvature.b(2,ak2)) + eta * Variables.aux  - 2.0 * adotoa * alpha; // M.Zaldarriaga's Ph.D. Thesis P.63
                                                                                    // Ma Bertschinger Eq(21d)
      //printf("\n\n\n %e %e %e %e %e %e",dgshear,ak2,Curvature.b(2,ak2),eta,Variables.aux,adotoa);

      //printf("\n\n%e %e %e %e %e %e %e",dgshear,ak2, Curvature.b(2,ak2),eta,Variables.aux,adotoa,alpha);
      //exit(1);
      
      //  Baryons.
      deltab = y[6];
      thetab = y[7];
      thetabdot = yprime[7];
      //  Photons.

      deltag = y[8];
      thetagdot = yprime[9];
      sheargdot = yprime[10] / 2.0;

      // Polarization term.
      polter = y[10] + y[9 + Variables.lmaxg] + y[11 + Variables.lmaxg];

      //printf("\n %e %e %e %e",Curvature.b(2, ak2), thetagdot, ak2, alphadot);
      //exit(1);
      coupldot = 8.0 * (Curvature.b(2, ak2) * thetagdot + Curvature.b(2, ak2) *  ak2 * alphadot) / 15.0;
      coupldot = coupldot - ak * 0.6 * Curvature.b(3, ak2) * (yprime[11] + yprime[10 + Variables.lmaxg] );
      coupldot = coupldot + ak * (1.0 - 0.4* Curvature.b(2, ak2)) * yprime[12 + Variables.lmaxg];

      polterdot = yprime[10] + yprime[9 + Variables.lmaxg] + yprime[11 + Variables.lmaxg];
      polterddot = coupldot - 0.3 * (GlobalArray.dopac[j] * polter + GlobalArray.opac[j] * polterdot);

      // Massless neutrinos.
      shearrdot = yprime[12 + 2 * Variables.lmaxg] / 2.0;

      //printf("\n\nSSS: %e %e %e",shearrdot,coupldot,polterddot);
      //exit(1);

      // Second derivative of expansion rate
      if (Variables.amnu == 0.0) // No massive neutrino
      {
         rhonudot = 0.0;
         shearnudot = 0.0;
      }
      else // Massive neutrino present
      {
         Neutrino.nuder(a, rhonu, ysanpass, adotoa, y + Variables.iq2 - 1, yprime + Variables.iq2 - 1);
         rhonudot = ysanpass[0];
         shearnudot = ysanpass[1];
      }

      // Dark energy
      omegavdyn = Variables.omegav * Darkenergy.dynrho(a);

      if(Variables.VPCflag == 1)  
         omegavdyn = omegavdyn * vpc.Variable_Lambda(a);


      weos = Darkenergy.wdyn_func(a);

      // Rho dot
      grhodot = (-vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a - 2.0 * (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + vpc.Variable_grhonr(a) * Variables.annunr * rhonu) / a2) * adotoa;
      grhodot = grhodot + vpc.Variable_grhonr(a) * Variables.annunr * rhonudot / a2 - ((1.0 + 3.0 * weos) * vpc.Variable_grhom(a) * omegavdyn * a2) * adotoa;

      if (Variables.dimflag == 1)
      {
         grho = pow((sqrt(vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / a2 + vpc.Variable_grhom(a) * Variables.omegav * vpc.Variable_Lambda(a) * a2) + sqrt(Variables.omegav * vpc.Variable_Lambda(a) * vpc.Variable_grhom(a)) * a), 2);
         sgrhooa2 = sqrt(grho / a2);
         rgrho = sqrt(vpc.Variable_grhom(a) * Variables.omegav * vpc.Variable_Lambda(a) + vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) / (a2 * a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / (a2 * a2));
         grhodot = (sgrhooa2 / rgrho) * (vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) * (-3 / a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) * (-4 / a2) + vpc.Variable_grhor(a) * Variables.annunr * rhonudot / (a2 * adotoa)) + 2 * grho;
      }

      adotoadot = grhodot / (6 * adotoa); // (ad/a)(d/dt)(ad/a) = (1/2)(d/dt)(ad/a)^2

      //  Derivative of the shear
      dgsheardot = (4.0 / 3.0) * (vpc.Variable_grhog(a) * sheargdot + Variables.annur * vpc.Variable_grhor(a) * shearrdot) / a2;
      dgsheardot = dgsheardot - 2.0 * adotoa * dgshear + Variables.annunr * vpc.Variable_grhonr(a) * shearnudot / a2;

      // Calculation of the sources
      alphaddot = -3 * dgsheardot / (2.0 * ak2) + etadot * Variables.aux  - 2.0 * adotoadot * alpha - 2.0 * adotoa * alphadot;

      s1 = etadot * Variables.aux  + alphaddot; // M.Zaldarriaga's Ph.D. Thesis Eq(3.38)
      s2 = 2 * alphadot;       // M.Zaldarriaga's Ph.D. Thesis Eq(3.38)

      d1 = GlobalArray.expmmu[j] * s1 + GlobalArray.vis[j] * (0.25 * deltag + s2 + polter / 16 / Curvature.b(2, ak2) + thetabdot / ak2 + (3.0 / 16.0) * polterddot / ak2 / Curvature.b(2, ak2));
      d1 = d1 + GlobalArray.dvis[j] * (alpha + thetab / ak2 + 3.0 / 8.0 * polterdot / ak2/Curvature.b(2, ak2) ) + GlobalArray.ddvis[j] * (3.0 / 16.0) * polter / ak2/Curvature.b(2, ak2);

      
      if(Variables.printISWseparate == 1){  	
          if(tau < Variables.taurrecombend)
            {
		d1 = d1 + GlobalArray.expmmu[j]*s1;
	  	zsanpass[4] = 0.0;
            }
            else 
                zsanpass[4] = GlobalArray.expmmu[j]*s1;
        }
      
      if (x > 0.0)
         dp1 = GlobalArray.vis[j] * (3.0 / 16.0) * polter / Curvature.b(2, ak2); //M.Zaldarriaga's Ph.D. Thesis Eq(3.38) Eq(3.37)
      else   //(x * x)/
         dp1 = 0.0;

      // lensing visibility function; approximate epoch of recombination
      phi1 = eta - adotoa * alpha;
      chi = tau0 - tau;
      chir = tau0 - Variables.taurmax;
      if (chi < chir) // lensing convergence, expmmu is an approximation one should integrate
      {               // visibility function over r(chi)/r(tau) but the error is harmless
         glens = (chir - chi) * chi / chir;
         dk1 = glens * ak2 * phi1 * GlobalArray.expmmu[j];
      }
      else
         dk1 = 0.0;

      //printf("\n %e %e %e %e %e %e",chi,chir,d1,dp1,dk1,phi1);
      //exit(1);
      zsanpass[0] = d1;
      zsanpass[1] = dp1;
      zsanpass[2] = dk1;
      zsanpass[3] = phi1;
   }





   void foutput(int n, double y[], double yprime[], int j, double tau0, double tau, double zsanpass[])
   {
      double x, a, a2, eta, etadot, alpha, alphadot, deltab, thetab;
      double thetabdot, deltag, thetagdot, sheargdot, polter, coupldot;
      double polterddot, shearrdot, rhonudot, shearnudot, omegavdyn;
      double weos, grhodot, grho, sgrhooa2, rgrho, adotoadot;
      double dgsheardot, alphaddot, s1, s2, chi, chir, glens;
      double polterdot, ysanpass[2];
      double d1, dp1, dk1, phi1;

      darkenergy Darkenergy;
      neutrino Neutrino;
      VPC vpc;

      d1 = 0.0;
      dp1 = 0.0;
      dk1 = 0.0;

      x = Variables.ak * (tau0 - tau);
      a = y[1];
      a2 = a * a;

      eta = y[3];
      etadot = yprime[3];

      alpha = (hdot + 6 * etadot) / (2.0 * Variables.ak2);                          // M.Zaldarriaga's Ph.D. Thesis P.63
      alphadot = -3 * dgshear / (2.0 * Variables.ak2) + eta - 2.0 * adotoa * alpha; // M.Zaldarriaga's Ph.D. Thesis P.63
                                                                                    // Ma Bertschinger Eq(21d)

      //  Baryons.
      deltab = y[6];
      thetab = y[7];
      thetabdot = yprime[7];
      //  Photons.

      deltag = y[8];
      thetagdot = yprime[9];
      sheargdot = yprime[10] / 2.0;

      // Polarization term.
      polter = y[10] + y[9 + Variables.lmaxg] + y[11 + Variables.lmaxg];

      coupldot = 8.0 * (thetagdot + Variables.ak2 * alphadot) / 15.0;
      coupldot = coupldot - Variables.ak * 0.6 * (yprime[11] + yprime[10 + Variables.lmaxg] + yprime[12 + Variables.lmaxg]);

      polterdot = yprime[10] + yprime[9 + Variables.lmaxg] + yprime[11 + Variables.lmaxg];
      polterddot = coupldot - 0.3 * (GlobalArray.dopac[j] * polter + GlobalArray.opac[j] * polterdot);

      // Massless neutrinos.
      shearrdot = yprime[12 + 2 * Variables.lmaxg] / 2.0;

      // Second derivative of expansion rate
      if (Variables.amnu == 0.0) // No massive neutrino
      {
         rhonudot = 0.0;
         shearnudot = 0.0;
      }
      else // Massive neutrino present
      {
         Neutrino.nuder(a, rhonu, ysanpass, adotoa, y + Variables.iq2 - 1, yprime + Variables.iq2 - 1);
         rhonudot = ysanpass[0];
         shearnudot = ysanpass[1];
      }

      // Dark energy
      omegavdyn = Variables.omegav * Darkenergy.dynrho(a);
      if(Variables.VPCflag == 1)  
         omegavdyn = omegavdyn * vpc.Variable_Lambda(a);

      weos = Darkenergy.wdyn_func(a);

      // Rho dot
      grhodot = (-vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a - 2.0 * (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * Variables.annur + vpc.Variable_grhonr(a) * Variables.annunr * rhonu) / a2) * adotoa;
      grhodot = grhodot + vpc.Variable_grhonr(a) * Variables.annunr * rhonudot / a2 - ((1.0 + 3.0 * weos) * vpc.Variable_grhom(a) * omegavdyn * a2) * adotoa;

      if (Variables.dimflag == 1)
      {
         grho = pow((sqrt(vpc.Variable_grhom(a) * (Variables.omegac + Variables.omegab) / a + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / a2 + vpc.Variable_grhom(a) * Variables.omegav * vpc.Variable_Lambda(a) * a2) + sqrt(Variables.omegav * vpc.Variable_Lambda(a) * vpc.Variable_grhom(a)) * a), 2);
         sgrhooa2 = sqrt(grho / a2);
         rgrho = sqrt(vpc.Variable_grhom(a) * Variables.omegav * vpc.Variable_Lambda(a) + vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) / (a2 * a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) / (a2 * a2));
         grhodot = (sgrhooa2 / rgrho) * (vpc.Variable_grhom(a) * (Variables.omegab + Variables.omegac) * (-3 / a) + (vpc.Variable_grhog(a) + vpc.Variable_grhor(a) * (Variables.annur + Variables.annunr * rhonu)) * (-4 / a2) + vpc.Variable_grhor(a) * Variables.annunr * rhonudot / (a2 * adotoa)) + 2 * grho;
      }

      adotoadot = grhodot / (6 * adotoa); // (ad/a)(d/dt)(ad/a) = (1/2)(d/dt)(ad/a)^2

      //  Derivative of the shear
      dgsheardot = (4.0 / 3.0) * (vpc.Variable_grhog(a) * sheargdot + Variables.annur * vpc.Variable_grhor(a) * shearrdot) / a2;
      dgsheardot = dgsheardot - 2.0 * adotoa * dgshear + Variables.annunr * vpc.Variable_grhonr(a) * shearnudot / a2;

      // Calculation of the sources
      alphaddot = -3 * dgsheardot / (2.0 * Variables.ak2) + etadot - 2.0 * adotoadot * alpha - 2.0 * adotoa * alphadot;

      s1 = etadot + alphaddot; // M.Zaldarriaga's Ph.D. Thesis Eq(3.38)
      s2 = 2 * alphadot;       // M.Zaldarriaga's Ph.D. Thesis Eq(3.38)

      d1 = GlobalArray.expmmu[j] * s1 + GlobalArray.vis[j] * (0.25 * deltag + s2 + polter / 16 + thetabdot / Variables.ak2 + (3.0 / 16.0) * polterddot / Variables.ak2);
      d1 = d1 + GlobalArray.dvis[j] * (alpha + thetab / Variables.ak2 + 3.0 / 8.0 * polterdot / Variables.ak2) + GlobalArray.ddvis[j] * (3.0 / 16.0) * polter / Variables.ak2;

      
      if(Variables.printISWseparate == 1){  	
          if(tau < Variables.taurrecombend)
            {
		d1 = d1 + GlobalArray.expmmu[j]*s1;
	  	zsanpass[4] = 0.0;
            }
            else 
                zsanpass[4] = GlobalArray.expmmu[j]*s1;
        }
      
      if (x > 0.0)
         dp1 = GlobalArray.vis[j] * (3.0 / 16.0) * polter / (x * x); //M.Zaldarriaga's Ph.D. Thesis Eq(3.38) Eq(3.37)
      else
         dp1 = 0.0;

      // lensing visibility function; approximate epoch of recombination
      phi1 = eta - adotoa * alpha;
      chi = tau0 - tau;
      chir = tau0 - Variables.taurmax;
      if (chi < chir) // lensing convergence, expmmu is an approximation one should integrate
      {               // visibility function over r(chi)/r(tau) but the error is harmless
         glens = (chir - chi) * chi / chir;
         dk1 = glens * Variables.ak2 * phi1 * GlobalArray.expmmu[j];
      }
      else
         dk1 = 0.0;

      zsanpass[0] = d1;
      zsanpass[1] = dp1;
      zsanpass[2] = dk1;
      zsanpass[3] = phi1;
   }

   /*************************************************************************************************************/
   /*		This subroutine computes the power spectra for mode Variables.ak of the scalar perturbations.							 */
   /*************************************************************************************************************/

   void powersflat(double ak, int in, double *apower)
   {
      double win, aknlog;
      //printf("\nHere in Scalar 808  %d\n In = %d\n",Variables.kcutflag,in);
      if (Variables.kcutflag == 0)
         win = 1.0;
      else
      {
         win = 2.0 * pow((ak / Variables.aksplit), 4);
         win = 2.0 * exp(-win) / (1.0 + exp(-win));
         if (Variables.kcutflag == -1)
            win = 1.0 - win;
         }

      // Normalize so tilt does not change power at pivot point k=0.05/Mpc
      aknlog = log(ak / 0.002);
      *apower = exp((Variables.an[in] - 1.0 + .5 * Variables.alphans[in]* aknlog + Variables.dalphansdlnk[in] * aknlog * aknlog / 6.) * aknlog);
      *apower = (*apower) * win;
      }

   void powersopen(double ak,int in,double *apower)
   {
      // This subroutine computes the power spectra
      // for mode ak of the scalar perturbations in the open case. 
      // Now it is set to a power law in k which is slightly more
      // complicated for beta(=ak).
  
      // Normalize so that tilt leaves the smae power 
      // at k=0.05 Mpc

      double anorm=0.05;
      double anorm2=anorm*anorm;
      double akn=ak/anorm;
      *apower = pow((ak*ak*Variables.r*Variables.r-Variables.kcurv*4.0),2)/pow((-Variables.kcurv+ak*Variables.r*ak*Variables.r),3)*Variables.r*Variables.r*ak*ak*akn*exp((Variables.an[in]-2.0)*0.5*log((-1.0*Variables.kcurv/Variables.r/Variables.r+ak*ak)/anorm2));
      *apower = *apower*100.0;
      }

};

void Fderivs(int n, double x, double y[], double yprime[])
{
   scalar Scalar;
   //printf("Here -- %e",Variables.omegak);
   if(fabs(Variables.omegak)<0.0001)
   {
      //printf(" --\n");
      Scalar.fderivs(n, x, y, yprime);
      //printf("----------------------------------------------------------------\n");
      }
   else
   {
      Scalar.fderivsCurved(n, x, y, yprime);
      }
}

#endif