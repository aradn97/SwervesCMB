#ifndef _DOUBLEINFLATION_H_
#define _DOUBLEINFLATION_H_

#include <stdio.h>
#include <math.h>
#include "initval.h"
#include "numericx.h"

static double km2mpc = 3.24077649e-20;
static double km2second = 1000.0/C;

double powerssH = 60.0;
double powersR = 5;
double powerss0 = 50;
double powersh00 = 72;

double findtheta(double rat,double R,double thetamin,double thetamax)		// thetamin,thetamax : upper and lower bound of theta
{
	double f1,f2,fm;
	double theta;
	theta = (thetamin+thetamax)/2;

	f1 = rat*pow(cos(thetamin),2*R*R/(R*R-1)) - pow(sin(thetamin),R*R/(R*R-1));
	f2 = rat*pow(cos(thetamax),2*R*R/(R*R-1)) - pow(sin(thetamax),R*R/(R*R-1));
	fm = rat*pow(cos(theta),2*R*R/(R*R-1)) - pow(sin(theta),R*R/(R*R-1));

	if(fabs(thetamax-thetamin)>0.0001)
	{
		if(f1*fm>0.0)
			theta = findtheta(rat,R,theta,thetamax);
		else
			theta = findtheta(rat,R,thetamin,theta);
		}
	return(theta);
	}

double generatephiS(double h0,double sH,	double R,double s0,double k)
{
	double kH;
	double sk;
	double thetak;
	double Hk; 
	double ml;
	double pi = 3.14159265;
	double phiSig;

	h0 = h0*km2second;										//Convert h0 in 1/Mpc
	kH = 2*pi*h0;
	sk = sH - log(k/kH);
	thetak = findtheta(sk/s0,R,0.0,pi/2);
	ml = 1.0; 													// This can be taken just as a normalization parameter

	Hk = (2.0/3.0)*ml*ml*sk*(1+(R*R-1)*sin(thetak)*sin(thetak));
	Hk = sqrt(Hk);	

	phiSig = 8*G/(9*pi)*Hk*Hk*sk;
	return(phiSig);
	}

double generatephiE(double h0,double sH,	double R,double s0,double k)
{
	double kH;
	double sk;
	double thetak;
	double Hk; 
	double ml;
	double pi = 3.14159265;
	double phiSig;

	h0 = h0*km2second;										//Convert h0 in 1/Mpc
	kH = 2*pi*h0;
	sk = sH - log(k/kH);
	thetak = findtheta(sk/s0,R,0.0,pi/2);
	ml = 1.0; 													// This can be taken just as a normalization parameter
	Hk = (2.0/3.0)*ml*ml*sk*(1+(R*R-1)*sin(thetak)*sin(thetak));
	Hk = sqrt(Hk);	

	phiSig = 2*G/(pi)*Hk*Hk/sk*(R*R*R*R/(cos(thetak)*cos(thetak)) + 1.0/(sin(thetak)*sin(thetak)));
	return(phiSig);
	}

double generatephiC(double h0,double sH,	double R,double s0,double k)
{
	double kH;
	double sk;
	double thetak;
	double Hk; 
	double ml;
	double pi = 3.14159265;
	double phiSig;

	h0 = h0*km2second;										//Convert h0 in 1/Mpc
	kH = 2*pi*h0;
	sk = sH - log(k/kH);

	thetak = findtheta(sk/s0,R,0.0,pi/2);
	ml = 1.0; 													// This can be taken just as a normalization parameter
	Hk = (2.0/3.0)*ml*ml*sk*(1+(R*R-1)*sin(thetak)*sin(thetak));
	Hk = sqrt(Hk);	
	
	phiSig = 4*G/(3*pi)*Hk*Hk*(R*R - 1);
	return(phiSig);
	}



double poweradiabatic(double k)
{
	double h00,sH;
	double R;
	double s0;

	sH = powerssH;
	R = powersR;
	s0 = powerss0;
	h00 = powersh00;
	
   return(generatephiS(h00,sH,R,s0,k));	
	}

double powerisotherm(double k)
{
	double h00,sH;
	double R;
	double s0;

	sH = powerssH;
	R = powersR;
	s0 = powerss0;
	h00 = powersh00;

	return(generatephiE(h00,sH,R,s0,k));	
	}

double powercross(double k)
{
	double h00,sH;
	double R;
	double s0;

	sH = powerssH;
	R = powersR;
	s0 = powerss0;
	h00 = powersh00;

	return(generatephiC(h00,sH,R,s0,k));	
	}









/****************************************************************************************************************************************/
//  	Reading power from the file 
/****************************************************************************************************************************************/

#define MAXTWOFIELDPOWER 10001

double twofieldpowersK[MAXTWOFIELDPOWER];
double twofieldpowersA[MAXTWOFIELDPOWER];
double twofieldpowersI[MAXTWOFIELDPOWER];
double twofieldpowersC[MAXTWOFIELDPOWER];

double twofieldpowersApr[MAXTWOFIELDPOWER];
double twofieldpowersIpr[MAXTWOFIELDPOWER];
double twofieldpowersCpr[MAXTWOFIELDPOWER];

int powerlow=0,powerup=MAXTWOFIELDPOWER;
int nmaxpow;

double powerA,powerI,powerC;


void readpowers()
{
	int i=0;
	char npowers[20];
	FILE *fppowers;
	
	printf("Enter the filename containing the initial power spectrum : \n");
	scanf("%s",npowers);
	printf("%s\n",npowers);
	
	fppowers = fopen(npowers,"r");
	printf("Hi");
	if(fppowers == NULL) 
	{
		printf("\nError : File don't exist / Problem in reading the file");
		exit(1);
		}

	while(!feof(fppowers))
		fscanf(fppowers,"%lf %lf %lf %lf",&twofieldpowersK[i],&twofieldpowersA[i],&twofieldpowersI[i],&twofieldpowersC[i++]);

	i--;
	nmaxpow = i;

	double d0hi = 1.0e40;
	double d0lo = 1.0e40;

	spline(twofieldpowersK,twofieldpowersA,i,d0lo,d0hi,twofieldpowersApr);
	spline(twofieldpowersK,twofieldpowersI,i,d0lo,d0hi,twofieldpowersIpr);
	spline(twofieldpowersK,twofieldpowersC,i,d0lo,d0hi,twofieldpowersCpr);

	printf("Number of lines : %d",i);
	
	}


void interpoletepowers(double ak)
{
	double h,a,b;
	int powermid;
	
	if(powerup == MAXTWOFIELDPOWER)
		powerup = nmaxpow;

	if(ak < twofieldpowersK[powerlow])
	{
		powerlow = 0;
		if(ak < twofieldpowersK[powerlow])
		{
			printf("The value of k is within the given table : \n Lowest k should be %e Mpc^-1\nSorry not possible to interpolate",ak);
			exit(1);
			//Check what can be returned
			}
		}

	if(ak > twofieldpowersK[powerup])
	{
		powerup = nmaxpow;
		if(ak > twofieldpowersK[powerup])
		{
			printf("The value of k is within the given table : \n Highest k should be %e Mpc^-1\nSorry not possible to interpolate",ak);
			exit(1);
			// Check How to tackle this .. 
			}
		}

	// Find the exact position in the table 

	while(powerup - powerlow > 1)
	{
		powermid = (powerup + powerlow)/2;
		if(ak > twofieldpowersK[powermid])
			powerlow = powermid;
		else 
			powerup=powermid;
		}

	// Interpolate
	
	//printf("\n%d %d %d",powerup,powerlow,powermid);
	h = twofieldpowersK[powerup] - twofieldpowersK[powerlow];

	if(h == 0.0)
	{
		printf("\nBad dataset for interpolation");
		exit(1);
		}

	a = (twofieldpowersK[powerup]-ak)/h;
	b = (ak - twofieldpowersK[powerlow])/h;	

	powerA=a*twofieldpowersA[powerlow]+b*twofieldpowersA[powerup]+((a*a*a-a)*twofieldpowersApr[powerlow]+(b*b*b-b)*twofieldpowersApr[powerup])*(h*h)/6.0;
	powerI=a*twofieldpowersI[powerlow]+b*twofieldpowersI[powerup]+((a*a*a-a)*twofieldpowersIpr[powerlow]+(b*b*b-b)*twofieldpowersIpr[powerup])*(h*h)/6.0;
	powerC=a*twofieldpowersC[powerlow]+b*twofieldpowersC[powerup]+((a*a*a-a)*twofieldpowersCpr[powerlow]+(b*b*b-b)*twofieldpowersCpr[powerup])*(h*h)/6.0;
	}




#endif