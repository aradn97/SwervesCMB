#ifndef CURVATURE_H_
#define CURVATURE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct evalujlpass{
   double chi,y1,y2;
} evalujlPass;

class curvature
{

double pi=3.141592653590;

double k;
public:
double curv;
numericx Numericx;

curvature(double curvature){
   curv = curvature;
}

double coshK(double x)
{
    if(k > 0.0) 
       return cos(x);
    else if (k < 0.0)
       return cosh(x);
    else
       return 1.0;
}

double sinhK(double x)
{
   //printf("\n\n%e\n",curv);
    if(curv > 0.0) 
       return sin(x);
    else if (curv < 0.0)
       return sinh(x);
    else
       return x;
}

double b(int l, double ak2)
{
   double b1 = 1.0-Variables.curv*(l*l-1.0)/ak2;

   if(fabs(Variables.omegak)<0.0001)
      return 1;

	if(b1 < 0.0) 
      return 0.0;
	return sqrt(b1);
}

//For Tensors

double bt1(int l, double betar)
{
   double bt1 = (1.0/l/l + 1.0/betar/betar)*(l*l-4.0);

   if(fabs(Variables.omegak)<0.0001)
      return 1;

	if (bt1 < 0.0) 
      return 0.0;
	return sqrt(bt1)*l;
}

double bt2(int l, double betar)
{
	double bt2 = (1.0/l/l + 1.0/betar/betar);

   if(fabs(Variables.omegak)<0.0001)
      return 1;
	
   if (bt2 < 0.0) 
      return 0.0;
	return sqrt(bt2)*(l*l-4.0);
}

double mod(double a,double p)
{
   return a - (int)(a/p)*p;
}


void symsource(double *s, double *s1, int iflip,int nstep)
{
   double stemp1,stemp2;
   for(int i=1; i<=iflip-1; i++)
   {
      stemp1=s[2*iflip-i];
      stemp2=s[i];
      s[2*iflip-i]=stemp1+stemp2;
      s1[2*iflip-i]=stemp1-stemp2;
      s[i]=s[2*iflip-i];
      s1[i]=-s1[2*iflip-i];
      }
   s[iflip]=2.0*s[iflip];
   s1[iflip]=0.0;
   for(int i=2*iflip; i<=nstep; i++)
         s1[i]=s[i];
   }

   void initchi(double beta,int j,double ll,double *y1,double *y2)
   {
      // For each a given beta and ll=l(j) find the initial
      // chi0 and give y1 and y2 to start the integration
      // of the bessel functions.

      double aux1, aux2;
      double a0,b0,h0,ho;
      int ibh, ibl, ib1;
      double betam1=1.0/beta;
      // Interpolate initial conditions

      //printf("\n%d",Variables.kcurv);
      //fflush(stdout);

      //Open models
      if(Variables.kcurv == -1) 
      {
         //printf("\nHIIIIIIIIIIIIIIi");
         //fflush(stdout);
         //printf("\n%e %e %e\n",beta,GlobalArray.abeta[Variables.ntot],GlobalArray.abeta[1]);
         //fflush(stdout);

         if((beta < GlobalArray.abeta[Variables.ntot]) && (beta > GlobalArray.abeta[1]))
         {
            // Beta is inside array
            ibl=(int)(Variables.dlnbeta*log(beta/Variables.betamin));
            ibh=ibl+1;
            ho=GlobalArray.abeta[ibh]-GlobalArray.abeta[ibl];
            a0=(GlobalArray.abeta[ibh]-beta)/ho;
            b0=(beta-GlobalArray.abeta[ibl])/ho;
            }
         else
         {
            // Extrapolate with a constant
            if (beta > GlobalArray.abeta[Variables.ntot])
            {
               ibl = Variables.ntot;
               ibh = Variables.ntot;
               a0=1.0;
               b0=0.0;
               ho=0.0;
               }
            // Extrapolate linearly
            else
            {
               ibl=1;
               ibh=ibl+1;
               ho=GlobalArray.abeta[ibh]-GlobalArray.abeta[ibl];
               a0=(GlobalArray.abeta[ibh]-beta)/ho;
               b0=(beta-GlobalArray.abeta[ibl])/ho;
               ho=0.0;
               }
            }
         aux2=a0*GlobalArray.ax[ibl][j] + b0*GlobalArray.ax[ibh][j] + ((a0*a0*a0-a0)*GlobalArray.axpr[ibl][j] +(b0*b0*b0-b0)*GlobalArray.axpr[ibh][j])*ho*ho/6.0; 
         aux1=exp(aux2)*betam1;
         Variables.chi0=log(aux1+sqrt(aux1*aux1 + 1.0));
         *y1=1.0e-6*aux1/ll;
         *y2=a0*GlobalArray.ay2[ibl][j]+b0*GlobalArray.ay2[ibh][j]+((a0*a0*a0-a0)*GlobalArray.ay2pr[ibl][j] +(b0*b0*b0-b0)*GlobalArray.ay2pr[ibh][j] )*ho*ho/6.0;
         }
      else
      {
         if(beta < GlobalArray.abetac[Variables.ntc[j]-Variables.Npoint+1][j])
         {
            // Beta is inside first part of array
            ibh=ll+1+Variables.ntc[j]-int(beta);
            ibl=ibh-1;
            ho=GlobalArray.abetac[ibh][j]-GlobalArray.abetac[ibl][j];
            a0=(GlobalArray.abetac[ibh][j]-beta)/ho;
            b0=(beta-GlobalArray.abetac[ibl][j])/ho;
            }
         else
         {
            if(beta < GlobalArray.abetac[1][j])
            {
               // Beta is inside second part of array
               ib1=int(Variables.dlnbeta*log(beta/Variables.betamin));
               ibh=Variables.ntot-ib1;
               ibl=ibh-1;
               if (beta > GlobalArray.abetac[ibl][j])
               {
                  ibl=ibl-1;
                  ibh=ibh-1;
                  }
               ho=GlobalArray.abetac[ibh][j]-GlobalArray.abetac[ibl][j];
               a0=(GlobalArray.abetac[ibh][j]-beta)/ho;
               b0=(beta-GlobalArray.abetac[ibl][j])/ho;
               }
            else
            {
               ibh=1;
               ibl=1;
               a0=1.0;
               b0=0.0;
               ho=0.0;
               }
            }

         aux2=a0*GlobalArray.axc[ibl][j]+b0*GlobalArray.axc[ibh][j]+((a0*a0*a0-a0)*GlobalArray.axcpr[ibl][j] +(b0*b0*b0-b0)*GlobalArray.axcpr[ibh][j])*ho*ho/6.0;
         aux1=exp(aux2)*betam1;
         Variables.chi0=asin(aux1);
         *y1=1.0e-6*aux1/ll;
         *y2  =a0*GlobalArray.ay2c[ibl][j]+b0*GlobalArray.ay2c[ibh][j]+((a0*a0*a0-a0)*GlobalArray.ay2cpr[ibl][j]+(b0*b0*b0-b0)*GlobalArray.ay2cpr[ibh][j])*ho*ho/6.0;
         }

      //printf("\nJIJIJIJIJIJIJIJIJ");
      //exit(1);         
      }


void evalujl(double *chi,double *y1,double *y2,double delchi,double beta2,double ap1)
{
   // Given y1 and y2, at chi, the function and its
   // derivative, move forward one step to chi+delchi.
   // At output chi, y1, y2 will be updated. 
   // K=-1 for open and +1 for closed. ap1=l*(l+1)
   // and beta2 is beta**2

   //evalujlPass temp;
   double dydchi1 = *y2;
   double dydchi2 = (ap1/pow(sinhK(*chi),2) - beta2)* (*y1);
   double hh = delchi*0.5;
   double h6 = delchi/6.0;
   double xh = *chi + hh;
   double yt1 = *y1 + hh*dydchi1;
   double yt2 = *y2 + hh*dydchi2;
   double dyt1 = yt2;
   double dyt2 = (ap1/pow(sinhK(xh),2) - beta2)*yt1;
   yt1 = *y1 + hh*dyt1;
   yt2 = *y2 + hh*dyt2;
   double dym1 = yt2;
   double dym2 = (ap1/pow(sinhK(xh),2) - beta2)*yt1;
   yt1 = *y1 + delchi*dym1;
   dym1 = dyt1 + dym1;
   yt2 = *y2 + delchi*dym2;
   dym2= dyt2 + dym2;
   dyt1= yt2;
   dyt2 = (ap1/pow(sinhK(*chi+delchi),2) - beta2)*yt1;
   *y1 = *y1 + h6*(dydchi1 + dyt1+2.0*dym1);
   *y2 = *y2 + h6*(dydchi2 + dyt2+2.0*dym2);      
   *chi = *chi + delchi;
   //return temp;
   }




//  For each a given beta and ll=l(j) find the initial
//  chi0 and give y1 and y2 to start the integration
//  of the bessel functions.


int tau2n(double tau,int nr)
{
   //  Given a tau finds the position in the atau0 array
   double an;
   int n;
   
   // printf("\nSKJKJKJKJ: %d",nr);fflush(stdout);
   for(int i=1; i<=nr; i++) 
   {
      // printf("\n%d %d",i,nr);fflush(stdout);
      // printf("\nFFFFF: %d %e %e",i, tau,GlobalArray.atau0[GlobalArray.nreg[i+1]]); fflush(stdout);

      if ((tau < GlobalArray.atau0[GlobalArray.nreg[i+1]]) && (tau >= GlobalArray.atau0[GlobalArray.nreg[i]]))
      {
         an=(tau-GlobalArray.atau0[GlobalArray.nreg[i]])/GlobalArray.dtaureg[i];
         n=GlobalArray.nreg[i]+(int)(an);
         return n;
         }
      }
   if (tau == GlobalArray.atau0[GlobalArray.nreg[nr+1]])
      n=GlobalArray.nreg[nr+1];
   else
   {
      printf("\nChecking Tau %e %e %d %d",tau,GlobalArray.atau0[GlobalArray.nreg[nr+1]],GlobalArray.nreg[nr+1],nr+1);
      printf("\ntau out of range in tau2n");
      printf("\n%e %e",tau,GlobalArray.atau0[GlobalArray.nreg[nr+1]]);
      exit(1);
      }
   return n;
   }


evalujlPass initujl(int ll,double ak,double tau,double *y1,double *y2,double tau0)
{
   // The bessel function is different from zero starting
   // at chi0, but we have to start the integration at recombination,
   // given by tau. So this routine advances from chi0,y1,y2 to
   // tau and also updates the values of y1 and y2.
   // If we eventually stop integrating the differential 
   // equation for the ujls and just call a subroutine of the form
   // ujl(chi,beta,l,K) the need to call initujl will 
   // disapear.

   //  TESTING: have a flag that sets if that particular ak is outputed

   double delchi;
   double beta=ak*Variables.r;
   double beta2=beta*beta;
   double betam1=1.0/beta;
   double chi;
   evalujlPass pass;

   // Integrating the diff. equation to get to tau
   double chi1=(tau0-tau)/Variables.r;

   //printf("\nCHI: %e %e",chi1,Variables.chi0);fflush(stdout);
   //   Finding the step in the for the ujl integration.
   if (chi1 > Variables.chi0)
   {
      if ((ll > 1) && (ll < 100))
         delchi=0.15*betam1;
      else
      {
         if(ll < 900)
            delchi=0.35*betam1;
         else
            delchi=0.50*betam1;
         }
         
      double sh=sinhK(Variables.chi0);
      double ch=coshK(Variables.chi0);
      double ap1 = 1.0*(ll*(ll+1));

      //printf("\nMU: %e %e %e",sh,ch,ap1);fflush(stdout);

      delchi=Numericx.max(delchi,0.30*sh/ch);

      //printf("\nMT: %e %e %e",sh,ch,ap1);fflush(stdout);

      int nstep=(int)((chi1-Variables.chi0)/delchi)+1;
      //printf("\nMS: %e %e %e",sh,ch,ap1);fflush(stdout);

      delchi=(chi1-Variables.chi0)/nstep;
      //printf("\nMR: %e %e %e",sh,ch,Variables.chi0);fflush(stdout);

      chi=Variables.chi0;

      //printf("\nMV: %e %e %e", sh, Variables.chi0, chi);fflush(stdout);
      for(int i=1; i <= nstep; i++)
      {                                        //  One step in the ujl integration
         //pass = 
         evalujl(&chi,y1,y2,delchi,beta2,ap1);   //  TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
//         chi = pass.chi;
//         *y1 = pass.y1;
//         *y2 = pass.y2;
         }
      //printf("\nMW: %e %e %e",sh,ch,ap1);fflush(stdout);         
      }
   pass.chi = chi;
   pass.y1 = *y1;
   pass.y2 = *y2;

   return pass;   
   }


int intopen1(int nstart, int nend,double dtau,double ak,int l,double tau0, double s2[], double sp2[], double sk2[],double ds2[],double dsp2[],double dsk2[],double *y1,double *y2,double *out1,double *out2,double *out3)
{
   // This subroutine integrates the source*ujl. 
   // It calculates ujl by integrating a second order
   // differential equation from initial values for calculating ujl.
   // atau0 is the array with the time where the sources are stored.
   // nstart and nend are the starting and finishing values of the
   // integration.
   // dtau is the spacing of the timestaps (they must be equally spaced)
   // s2, sp2 and sk2 are the sources. ds2, dp2, dsk2 are their second 
   // derivatives for the spline interpolation. 
   // ak is the wavevector, l is the multipole.
   // tau0 is the time today.
   // r is the curvature radius.

   // dxsource is the spacing in x of the sampling of the source.
   // dxmax is the maximun step for which the ujl integration converges
   // correctly.


   evalujlPass pass;
   numericx Numericx;   
   double aux1, aux2, beta, betam1, betam2, beta2, ap1, chi;
   double sh, shml, ujl, shm1;
   double dchimax, dchisource, delchi, dtau2;
   double st,sp,sk;   
   double Deltachi, taui;
   int is;
   double xa,xb;
   
   //printf("\n-----(N1)-"); fflush(stdout);
   //printf("\n ---------------- %d %d",nstart,nend); fflush(stdout);
   if(nstart == nend)
   {
      *out1=0.0;
      *out2=0.0;
      *out3=0.0;
      return 0;
      }

   //printf("\n-----(N2)-");fflush(stdout);

      aux1=1.0*Variables.r/dtau;
      aux2=(tau0-GlobalArray.atau0[nstart])/dtau + nstart;

      beta=ak*Variables.r;
      betam1=1.0/beta;
      beta2=beta*beta;
      betam2=betam1*betam1;
      ap1=1.0*l*(l+1);

      chi=(tau0-GlobalArray.atau0[nend])/Variables.r;
      sh=sinhK(chi);
      shm1=1.0/sh;
      ujl=*y1*shm1;

      //printf("\nTHIS:%e %e %e %e",chi,sh,shm1,*y1);fflush(stdout);
      //printf("\n-----(N3)-");
      
      if ((l > 1) && (l < 100))
         dchimax=0.15*betam1;
      else
      {
         if (l < 900) 
            dchimax=0.35*betam1;
         else
            dchimax=0.50*betam1;
         }

      // printf("\n-----(N4)-");fflush(stdout);
      dchimax=Numericx.min(dchimax,fabs(0.30*sh/(coshK(chi)+0.10)));
      dchisource=dtau/Variables.r;

      // printf("\n-----(N4a)-");fflush(stdout);
      // printf("\n-----(N4b):%d %e %e",nend,ujl,s2[nend]);fflush(stdout);
      // printf("\n-----(N4c):%e %e %e",betam2,shm1,shm1);fflush(stdout);

      *out1=0.50*ujl*s2[nend];
      *out2=0.50*ujl*sp2[nend]*betam2*shm1*shm1;
      *out3=0.50*ujl*sk2[nend];

      // printf("\n-----(N5axx)-");fflush(stdout);      
      // printf("\n-----(N5bxx)-");fflush(stdout);            
      // printf("\n-----(N5)-");fflush(stdout);

      //printf("\n%d %e %e %e %e %e",nend,ujl,s2[nend],*out1,*out2,*out3);

      if (dchisource > dchimax)
      {
         // Interpolate the source.
         delchi=dchimax;
         Deltachi=(GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/Variables.r;
         int nstep=(int)(Deltachi/delchi)+1;
         delchi=Deltachi/nstep;
         dtau2=Variables.r*Variables.r*delchi*delchi;
         //printf("\nIF:");
         //printf("\n-----(N3-if)-");fflush(stdout);
         for(int i=1;i<=nstep;i++)         // One step in the ujl integration
         {
            //pass = 
            evalujl(&chi,y1,y2,delchi,beta2,ap1);
            //chi = pass.chi;
            //y1 = pass.y1;
            //y2 = pass.y2;

            shm1=1.00/sinhK(chi);
            ujl=*y1 * shm1;


            // Interpolate the source
            taui=aux2-aux1*chi;
            is=(int)taui;

            //printf("\n%d %e %e %e %e",is,taui,aux2,aux1,chi);
            xb=taui-is;
            xa=1.0-xb;
            st=xa*s2[is]+xb*s2[is+1]+((xa*xa*xa-xa)*ds2[is]+(xb*xb*xb-xb)*ds2[is+1])*dtau2/6.0;
            sp=xa*sp2[is]+xb*sp2[is+1]+((xa*xa*xa-xa)*dsp2[is]+(xb*xb*xb-xb)*dsp2[is+1])*dtau2/6.0;
            sk=xa*sk2[is]+xb*sk2[is+1]+((xa*xa*xa-xa)*dsk2[is]+(xb*xb*xb-xb)*dsk2[is+1])*dtau2/6.0;

            *out1 = *out1+st*ujl;
            *out2 = *out2+sp*ujl*betam2*shm1*shm1;
            *out3 = *out3+sk*ujl;
            }

         *out1 = *out1-0.50*st*ujl;
         *out2 = *out2-0.50*sp*ujl*betam2*shm1*shm1;
         *out3 = *out3-0.50*sk*ujl;

         *out1 = *out1*delchi*Variables.r;
         *out2 = *out2*delchi*Variables.r;
         *out3 = *out3*delchi*Variables.r;
         }
      else
      {
         //printf("\nELSE:");
         // Compute jl at values where the source is stored. 
         delchi=dchisource;


         for(int i=(nend-1); i>=nstart; i--)
         {
            // One step in the ujl integration
            //pass = 
            //printf("\n%d %e %e",i,s2[i],sp2[i]);
            evalujl(&chi,y1,y2,delchi,beta2,ap1);
            //chi = pass.chi;
            //y1 = pass.y1;
            //y2 = pass.y2;

            shm1=1.00/sinhK(chi);
            ujl=*y1*shm1;

            *out1 = *out1+s2[i]*ujl;
            *out2 = *out2+sp2[i]*ujl*betam2*shm1*shm1;
            *out3 = *out3+sk2[i]*ujl;
            //printf("\nYU %d %d %e %e %e %e %e %e %e",l,i,ak,tau0,s2[i],ujl,*y1,shm1,*out1);            
            }
         //exit(1);   
         *out1 = *out1-0.50*s2[nstart]*ujl;
         *out2 = *out2-0.50*sp2[nstart]*ujl*betam2*shm1*shm1;
         *out3 = *out3-0.50*sk2[nstart]*ujl;

         *out1 = *out1 * dtau;
         *out2 = *out2 * dtau;
         *out3 = *out3 * dtau;

         //printf("\n%d %e %e %e %e %e",nstart,ujl,s2[nstart],*out1,*out2,*out3);
      }

      //printf("\nIntopen1\n%e %e %e",*out1,*out2,*out3);
      return 0;
   }

int intopen2(int nstart, int nend,double dtau,double ak,int l,double tau0, double s2[], double sp2[], double sk2[],double ds2[],double dsp2[],double dsk2[],double *y1,double *y2,double *out1,double *out2,double *out3)
{
   // This subroutine integrates the source*ujl. 
   // It approximates ujl by an asintotic formula using y1 and y2
   // to match continuously with the exact ujl. 
   // atau0 is the array with the time where the sources are stored.
   // nstart and nend are the starting and finishing values of the
   // integration.
   // dtau is the spacing of the timestaps (they must be equally spaced)
   // s2 and sp2 are the sources. ds2 dp2 are their second derivatives for
   // the spline interpolation. 
   // ak is the wavevector, l is the multipole.
   // tau0 is the time today.
   // r is the curvature radius.


      double aux1, aux2, beta, betam1, betam2, beta2;
      double dchimax, dchisource, aux01;
      double dtau2,g2,chi,shm1,dtaup;
      double amp,phi,delta,delchi,taui;
      double ujl,delphi,xa,xb;
      double st,sp,sk;
      int is,nstp;

      numericx Numericx;

      if (nstart == nend)
      {
         *out1=0.0;
         *out2=0.0;
         *out3=0.0;
         return 0;
         }

      aux1=1.0*Variables.r/dtau;
      aux2=(tau0-GlobalArray.atau0[nstart])/dtau + nstart;
      dtau2=dtau*dtau;

      betam1=1.00/ak/Variables.r;
      betam2=betam1*betam1;
      g2=(l*(l+1))/pow((ak*Variables.r),2);

      chi=(tau0-GlobalArray.atau0[nend])/Variables.r;
      shm1=1.00/sinhK(chi);
      amp=betam1/pow((1.00-g2*shm1*shm1),0.25);
      
      aux01 = *y1/amp;
      if (fabs(aux01) >= 1.)
      {
         *y1 = (*y1) * fabs(amp/(*y1));
         aux01 = (*y1)/fabs(*y1);
         }

      if (*y2 >= 00) 
         delta=asin(aux01);
      else
         delta=pi-asin(aux01);


      dchimax=0.50*betam1;
      dchisource=dtau/Variables.r;

   // There is no need to integrate modes that oscillate
   // too much in the time the source remains constant, so return.

      if (dchimax < (dchisource*0.50))
      {
         *out1=0.00;
         *out2=0.00;
         *out3=0.00;
         }
      
      delchi=Numericx.min(dchimax,dchisource);
      dtaup=Variables.r*delchi;

      nstp=(int)((GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/dtaup)+1;
      dtaup=(GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/nstp;
      delchi=dtaup/Variables.r;

      ujl = *y1 * shm1;
      *out1 = 0.50*s2[nend]*ujl;
      *out2 = 0.50*sp2[nend]*ujl*betam2*shm1*shm1;
      *out3 = 0.50*sk2[nend]*ujl;

      delphi=ak*dtaup;

   // In models where the source has been flipped because
   // chi to recombination is more than pi/2 make sure ujl 
   // has the correct parity. This means that when we switch 
   // to WKB the ujl is discontineous. This does not cause significant
   // error. Instead if the parity is not correct and pi/2 falls inside
   // recombination we run into big trouble because the cancelation of the
   // velocity term does not occur (because of the integration by parts) so
   // modes with a significant velocity contribution are wrong. So we need 
   // to make sure that the parity is correct. In models where pi/2 is not 
   // inside recombiantion, all is irrelevant. 
 
      if(Variables.ichiflag == 1)
      {
         beta=ak*Variables.r;
         delta=beta*chi-l*pi/2.00;
         delta=mod(delta,2.0*pi);
         }

      phi=delta;

      for(int i=1;i<=nstp;i++)
      {
         chi=chi+delchi;
         shm1=1.00/sinhK(chi);
         phi=phi+delphi;
         ujl=betam1*shm1*sin(phi)/pow((1.00-g2*shm1*shm1),0.25);

         // Interpolate the source
         taui=aux2-aux1*chi;
         is=(int)taui;
         xb=taui-is;
         xa=1.0-xb;
         st=xa*s2[is]+xb*s2[is+1]+((xa*xa*xa-xa)*ds2[is]+(xb*xb*xb-xb)*ds2[is+1])*dtau2/6.0;
         sp=xa*sp2[is]+xb*sp2[is+1]+((xa*xa*xa-xa)*dsp2[is]+(xb*xb*xb-xb)*dsp2[is+1])*dtau2/6.0;
         sk=xa*sk2[is]+xb*sk2[is+1]+((xa*xa*xa-xa)*dsk2[is]+(xb*xb*xb-xb)*dsk2[is+1])*dtau2/6.0;

         *out1 = *out1+ujl*st;
         *out2 = *out2+ujl*sp*betam2*shm1*shm1;
         *out3 = *out2+ujl*sk;
         }

      *out1 = *out1-0.50*ujl*s2[nstart];
      *out1 = *out1 * dtaup;
      *out2 = *out2-0.50*ujl*sp2[nstart]*betam2*shm1*shm1;
      *out2 = *out2 * dtaup;
      *out3 = *out3-0.50*ujl*sk2[nstart];
      *out3 = *out3 * dtaup;
      
      *y1=ujl*sinhK(chi);
      *y2=cos(phi);

      //iphase=1;
      return 0;
   }

double intopen1t(int nstart, int nend, double dtau,double ak,int l,double tau0, double s2[], double se2[], double sb2[],double ds2[],double dse2[],double dsb2[],double *y1,double *y2,double *out1,double *out2,double *out3)
{
   // This subroutine integrates the source*ujl. 
   // It calculates ujl by integrating a second order
   // differential equation from initial values for calculating ujl.
   // atau0 is the array with the time where the sources are stored.
   // nstart and nend are the starting and finishing values of the
   // integration.
   // dtau is the spacing of the timestaps (they must be equally spaced)
   // s2, se2 and sb2 are the sources. ds2 dp2 are their second derivatives for
   // the spline interpolation. 
   // ak is the wavevector, l is the multipole.
   // tau0 is the time today.
   // r is the curvature radius.

   //   integer nstart,nend

   // dxsource is the spacing in x of the sampling of the source.
   // dxmax is the maximun step for which the ujl integration converges
   // correctly.

      double aux1, aux2, beta, betam1, beta2;
      double xa, xb, sh, chi;
      double taui, dchimax, ujl;
      double deltachi,delchi,dchisource,ap1,dtau2;
      double st,sb,se;
      int is,nstep;
      evalujlPass pass;
      numericx Numericx;

      if (nstart == nend) 
      {
         *out1=0.0;
         *out2=0.0;
         *out3=0.0;
         return 0;
         }

      aux1=1.0*Variables.r/dtau;
      aux2=(tau0-GlobalArray.atau0[nstart])/dtau + nstart;

      beta=ak*Variables.r;
      betam1=1.0/beta;
      beta2=beta*beta;
      ap1=1.0*(l*(l+1));


      chi=(tau0-GlobalArray.atau0[nend])/Variables.r;
      sh=sinhK(chi);
      ujl=*y1/sh;

      if ((l > 1) && (l < 300))
      {
         if (l < 60)
            dchimax=0.15*betam1;
         else
         {
            dchimax=0.06*betam1;
            if (l >= 600) dchimax=0.03*betam1;
            }
         }
      else
      {
         if (l < 900) 
            dchimax=0.35*betam1;
         else
            dchimax=0.50*betam1;
         }

      dchimax = Numericx.min(dchimax,fabs(0.30*sh/(coshK(chi)+0.10)));
      dchisource=dtau/Variables.r;

      *out1 = 0.50*ujl*s2[nend];
      *out2 = 0.50*ujl*se2[nend];
      *out3 = 0.50*ujl*sb2[nend];

      if (dchisource > dchimax)
      {
         // Interpolate the source.
         delchi=dchimax;
         deltachi=(GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/Variables.r;
         nstep=(int)(deltachi/delchi)+1;
         delchi=deltachi/nstep;
         dtau2=Variables.r*Variables.r*delchi*delchi;

         for(int i=1; i<=nstep; i++)
         {
            // One step in the ujl integration
            //pass = 
            evalujl(&chi,y1,y2,delchi,beta2,ap1);
            //chi = pass.chi;
            //y1 = pass.y1;
            //y2 = pass.y2;

            ujl=*y1/sinhK(chi);

            // Interpolate the source
            taui=aux2-aux1*chi;
            is = (int)taui;
            xb=taui-is;
            xa=1.00-xb;
            st=xa*s2[is]+xb*s2[is+1]+((xa*xa*xa-xa)*ds2[is]+(xb*xb*xb-xb)*ds2[is+1])*dtau2/6.0;
            se=xa*se2[is]+xb*se2[is+1]+((xa*xa*xa-xa)*dse2[is]+(xb*xb*xb-xb)*dse2[is+1])*dtau2/6.0;
            sb=xa*sb2[is]+xb*sb2[is+1]+((xa*xa*xa-xa)*dsb2[is]+(xb*xb*xb-xb)*dsb2[is+1])*dtau2/6.0;

            *out1 = *out1+st*ujl;
            *out2 = *out2+se*ujl;
            *out3 = *out3+sb*ujl;
            }

         *out1 = *out1-0.50*st*ujl;
         *out2 = *out2-0.50*se*ujl;
         *out3 = *out3-0.50*sb*ujl;

         *out1 = *out1*delchi*Variables.r;
         *out2 = *out2*delchi*Variables.r;
         *out3 = *out3*delchi*Variables.r;
         }
      else
      {
         // Compute jl at values where the source is stored. 
         delchi=dchisource;
         for(int i=(nend-1); i>=nstart; i--)
         {
            // One step in the ujl integration
            //pass = 
            evalujl(&chi,y1,y2,delchi,beta2,ap1);
            //chi = pass.chi;
            //y1 = pass.y1;
            //y2 = pass.y2;

            ujl=*y1/sinhK(chi);
            *out1 = *out1+s2[i]*ujl;
            *out2 = *out2+se2[i]*ujl;
            *out3 = *out3+sb2[i]*ujl;
            }

         *out1 = *out1-0.50*s2[nstart]*ujl;
         *out2 = *out2-0.50*se2[nstart]*ujl;
         *out3 = *out3-0.50*sb2[nstart]*ujl;

         *out1 = *out1*dtau;
         *out2 = *out2*dtau;
         *out3 = *out3*dtau;
         }
   return 0;         
   }      

double intopen2t(int nstart, int nend, double dtau,double ak,int l,double tau0, double s2[], double se2[], double sb2[],double ds2[],double dse2[],double dsb2[],double *y1,double *y2,double *out1,double *out2,double *out3)
{
   // This subroutine integrates the source*ujl. 
   // It approximates ujl by an asintotic formula using y1 and y2
   // to match continuously with the exact ujl. 
   // atau0 is the array with the time where the sources are stored.
   // nstart and nend are the starting and finishing values of the
   // integration.
   // dtau is the spacing of the timestaps (they must be equally spaced)
   // s2, se2 and sb2 are the sources. ds2 dp2 are their second derivatives for
   // the spline interpolation. 
   // ak is the wavevector, l is the multipole.
   // tau0 is the time today.
   // r is the curvature radius.


      double aux1, aux2, beta, betam1, beta2;
      double taui,delta,phi;
      double xa,xb;
      double st,se,sb;
      double ujl,delphi, dtaup;
      double dchimax,dchisource,chi;
      double dtau2,g2,shm1,amp,aux01,delchi;
      int is, nstp;
      numericx Numericx;

      if (nstart == nend)
      {
         *out1=0.00;
         *out2=0.00;
         *out3=0.00;
         return 0;
         }

      aux1=1.0*Variables.r/dtau;
      aux2=(tau0-GlobalArray.atau0[nstart])/dtau + nstart;
      dtau2=dtau*dtau;

      betam1=1.0/ak/Variables.r;
      g2=(l*(l+1))/pow((ak*Variables.r),2);

      chi=(tau0-GlobalArray.atau0[nend])/Variables.r;
      shm1=1.00/sinhK(chi);
      amp=betam1/pow((1.00-g2*shm1*shm1),0.25);
      
      aux01=*y1/amp;
      if (fabs(aux01) >= 1.0)
      {
         *y1=*y1*fabs(amp/ *y1);
         aux01=*y1/fabs(*y1);
         }

      if (*y2 >= 00)
         delta=asin(aux01);
      else
         delta=pi-asin(aux01);

      
      dchimax=0.50*betam1;

      dchisource=dtau/Variables.r;

      delchi=Numericx.min(dchimax,dchisource);
      dtaup=Variables.r*delchi;

      nstp=(int)((GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/dtaup)+1;
      dtaup=(GlobalArray.atau0[nend]-GlobalArray.atau0[nstart])/nstp;
      delchi=dtaup/Variables.r;


      ujl=*y1*shm1;
      *out1=0.50*s2[nend]*ujl;
      *out2=0.50*se2[nend]*ujl;
      *out3=0.50*sb2[nend]*ujl;

      delphi=ak*dtaup;

      if (Variables.ichiflag == 1)
      {
         beta=ak*Variables.r;
         delta=beta*chi-l*pi/2.00;
         delta=mod(delta,2.00*pi);
         }

      phi=delta;

      for(int i=1;i<=nstp;i++)
      {
         chi=chi+delchi;
         shm1=1.00/sinhK(chi);
         phi=phi+delphi;
         ujl=betam1*shm1*sin(phi)/pow((1.00-g2*shm1*shm1),0.25);

         // Interpolate the source
         taui=aux2-aux1*chi;
         is=int(taui);
         xb=taui - is;
         xa=1.00-xb;
         st=xa*s2[is]+xb*s2[is+1]+((xa*xa*xa-xa)*ds2[is]+(xb*xb*xb-xb)*ds2[is+1])*dtau2/6.0;
         se=xa*se2[is]+xb*se2[is+1]+((xa*xa*xa-xa)*dse2[is]+(xb*xb*xb-xb)*dse2[is+1])*dtau2/6.0;
         sb=xa*sb2[is]+xb*sb2[is+1]+((xa*xa*xa-xa)*dsb2[is]+(xb*xb*xb-xb)*dsb2[is+1])*dtau2/6.0;

         *out1 = *out1+ujl*st;
         *out2 = *out2+ujl*se;
         *out3 = *out3+ujl*sb;
         }

      *out1 = *out1-0.50*ujl*s2[nstart];
      *out1 = *out1*dtaup;
      *out2 = *out2-0.50*ujl*se2[nstart];
      *out2 = *out2*dtaup;
      *out3 = *out3-0.50*ujl*sb2[nstart];
      *out3 = *out3*dtaup;
      
      *y1=ujl*sinhK(chi);
      *y2=cos(phi);

      return 0;
      }


void initujl0()
{
   // This subroutine reads the ujl files from disk and 
   // initializes other variables needed in CMBFAST.

   int lmofile,l0file,ntotc,ist;
   double betamaxfile,aist;
   double d0hi = 1.0e40, d0lo = 1.0e40;   
   FILE *fp1;
   
   fp1 = fopen("myujl.dat", "r"); 

   fscanf(fp1,"%d", &lmofile);
   fscanf(fp1,"%lf %lf %lf", &betamaxfile, &Variables.betamin, &Variables.dlnbeta);
   fscanf(fp1,"%d", &Variables.Npoint);

  
   if (Variables.lmo > lmofile)
   {
      printf("\nYou have entered a lmax");
      printf("\ninconsistent with those in the file");
      printf("\n%d",&lmofile);
      printf("\nYou will have to start again");
      exit(1);
      }

   // Checking if the lvalues.inc file used to build jl file
   // is the same as the one in the code.
   fscanf(fp1,"%d", &l0file);

   for(int j=1; j<=l0file; j++)
   {
      fscanf(fp1,"%d", &Variables.lfile[j]);
      }
    

   for(int j=1;j<=Variables.l0;j++)
   {
      if (Variables.l[j] != Variables.lfile[j])
      {
         printf("\nlvalues.inc file used to build jl file");
         printf("\nand the one in the code differ.");
         printf("\nYou must use the same one");
         printf("\n%d %d %d",j,Variables.l[j],Variables.lfile[j]);
         exit(1);
         }
      }
     //exit(1);
   // reading  uj_l information 
   // remember to create ujl.dat with ujlgen.f first using 
   // correct lmax 


   //printf("\n\n %e %e %e",betamaxfile,Variables.betamin,Variables.dlnbeta);
   // Number of betas for open models.
   Variables.ntot = (int)(log(betamaxfile/Variables.betamin)*Variables.dlnbeta);
   //printf("\n\n%d",Variables.ntot);exit(1);

   for(int i=1; i <= Variables.ntot; i++)
      GlobalArray.abeta[i]=Variables.betamin*exp((i-1)/Variables.dlnbeta);
      //printf("\n%e %e", Variables.betamin,Variables.dlnbeta);
      //}
      
//exit(1);
   // Read open models

   for(int j=1; j<=Variables.l0; j++)
   {
      for(int i=1; i<=Variables.ntot; i++)
      {
         fscanf(fp1,"%lf", &GlobalArray.ax[i][j]);
         fscanf(fp1,"%lf", &GlobalArray.ay2[i][j]);
         //printf("\n%d %d %e %e",j,i,GlobalArray.ax[i][j],GlobalArray.ay2[i][j]);
         }
      //fflush(stdout);exit(1);  
      // Closed models

      // The list of betas for closed models will be
      // l+1,l+2,...,l+Npoint,beta(ist+1),beta(ist+2),....,beta(ntot)
      // where ist is given by

      aist=1.0+Variables.dlnbeta*log((Variables.l[j] + Variables.Npoint+1)/Variables.betamin);
      ist=(int)aist;
      ntotc=Variables.ntot-ist + Variables.Npoint;

      // Note that ntotc depends on l. We keep this numbers in an array
      Variables.ntc[j]=ntotc;
//exit(1);
      for(int i=1; i<=ntotc; i++)
      {
         // So that the betas for this l are
//         printf("\n%d %d",i,ntotc); fflush(stdout);

         if (i <= (ntotc - Variables.Npoint))
         {
            GlobalArray.abetac[i][j]=GlobalArray.abeta[Variables.ntot-i+1];
            GlobalArray.abetac[i][j]=((int)(GlobalArray.abetac[i][j]));
            }
         else
            GlobalArray.abetac[i][j]=1.0*(Variables.l[j]+1+ntotc-i);

         fscanf(fp1,"%lf", &GlobalArray.axc[i][j]);
         fscanf(fp1,"%lf", &GlobalArray.ay2c[i][j]);
         //printf("\n%e %e",GlobalArray.axc[i][j],GlobalArray.ay2c[i][j]); fflush(stdout);
         }
      }
   
   //printf("\n%d %d",Variables.lmo, lmofile);   
   //exit(1);   
   fclose(fp1);
   
   //  get the interpolation matrices for starting value and deriv.

   //printf("\n%e %e %e",GlobalArray.abetac[0][0],GlobalArray.axc[0][0],GlobalArray.axcpr[0][0]);
   //exit(1);
   Numericx.spline2D(GlobalArray.abeta,GlobalArray.ax,Variables.ntot,Variables.l0,d0lo,d0hi,GlobalArray.axpr);
   Numericx.spline2D(GlobalArray.abeta,GlobalArray.ay2,Variables.ntot,Variables.l0,d0lo,d0hi,GlobalArray.ay2pr);
   Numericx.spline2DLineWise(GlobalArray.abetac,GlobalArray.axc,Variables.ntc,Variables.l0,d0lo,d0hi,GlobalArray.axcpr);
   //Numericx.spline2DLineWise(GlobalArray.abetac,GlobalArray.axc,Variables.ntc,Variables.l0,d0lo,d0hi,GlobalArray.axcpr);   
   //exit(1);
   Numericx.spline2DLineWise(GlobalArray.abetac,GlobalArray.ay2c,Variables.ntc,Variables.l0,d0lo,d0hi,GlobalArray.ay2cpr);
   //exit(1);
   }

};



#endif