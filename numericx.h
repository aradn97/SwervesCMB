#ifndef	_NUMERICX_H_
#define	_NUMERICX_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "initval.h"

/************************************************************************************************************************************/
/*		Min function will take two values and calculate the menimum of those two values. This function has been used several times in */
/*		our programe																						 																*/
/************************************************************************************************************************************/

double min(double compx, double compy)
{
	return (compx<compy)?compx:compy;
	}

int minint(int compx, int compy)
{
	return (compx<compy)?compx:compy;
	}

/************************************************************************************************************************************/
/*		Max function will take two values and calculate the maximum of those two values. This function has been used several times in */
/*		our programe																						 																*/
/************************************************************************************************************************************/

double max(double compx, double compy)
{
	return (compx>compy)?compx:compy;
	}

int maxint(int compex, int compey)	
{
	return (compex>compey)?compex:compey;
	}



/************************************************************************************************************************************/
//		Numerical Subroutine Variables																															   //
//    Developer's Guide : If "Error (splder):Array overflow" occured then increase the array size. Also make the necessery changes  //
//    in splini() and splder() 			          																												//
/************************************************************************************************************************************/
double gg[40002];	

// g is dummy variable. 																					 															//
// Used by : 1) splini() 2) splder() 																	 															//	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/************************************************************************************************************************************/
/*		splini() should be called everytime before calling the function splder()																		*/
/************************************************************************************************************************************/

void splini()
{
	gg[1]=0.0;
	for(int i=2;i<=40001;i++)
		gg[i]=1.0/(4.0-gg[i-1]);
	}


/************************************************************************************************************************************/
/*		Splder fits a cubic spline to y and returns the first derivatives at the grid points in dy.												*/
/*   	Dy is equivalent to a 4th-order Pade difference formula for dy/di.																				*/
/************************************************************************************************************************************/

void splder(double y[],double dy[],int n)
{
	double f[40002];
	int nl=n-1;
	
	if(nl>40000) printf("Error (splder):Array overflow");
	
	f[1]=(-10.0*y[1]+15.0*y[2]-6.0*y[3]+y[4])/6.0;
   f[n]=(10.0*y[n]-15.0*y[nl]+6.0*y[n-2]-y[n-3])/6.0;
   
   for(int i=2;i<=nl;i++)
     f[i]=gg[i]*(3.0*(y[i+1]-y[i-1])-f[i-1]);
	
	dy[n]=f[n];
	
	for(int i=nl;i>=1;i--)
		dy[i]=f[i]-gg[i]*dy[i+1];	   	
	}
	

/************************************************************************************************************************************/
/* 	Splint integrates a cubic spline, providing the ouput value																							*/
/*	 	z = integral from 1 to n of s(i)di, where s(i) is the spline fit to y(i).																		*/
/************************************************************************************************************************************/

double splint(double y[],int n)
{
	int nl=n-1;
	double dy1=0.0;
	double dyn=(11.0*y[n]-18.0*y[nl]+9.0*y[n-2]-2.0*y[n-3])/6.0;
	double z=0.5*(y[1]+y[n-1])+(dy1-dyn)/12.0;
	for(int i=2;i<nl+1;i++)
	{
		z=z+y[i];
		}
	return z;
	}	

/*******************************************************************************************************************/
//  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
//  Spline 
//  For details check this : http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f3-3.pdf
//
// Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi= f (xi), with x1 < x2 < . . . < xN , 
// and given values yp1 and ypn for the ?rst derivative of the interpolating function at points 1 and n, respectively, 
// this routine returns an array y2(1:n) of length n which contains the second derivatives of the interpolating function 
// at the tabulated points xi. If yp1 and/or ypn are equal to 1 × 10^30 or larger, the routine is signaled to set
// the corresponding boundary condition for a natural spline, with zero second derivative on that boundary.
//
// Parameter: NMAX is the largest anticipated value of n
/*******************************************************************************************************************/

void spline(double x[],double y[],int n,double yp1,double ypn,double y2[])
{
	double qn,un,p,sig;
	double u[100010];
   
   if(yp1 > .99e30)
   {
	   y2[1]=0.0;
      u[1]=0.0;
      }
	
	else
	{
		y2[1]=-0.50;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
      }
      
	for(int i=2;i<=n-1;i++)
	{
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		}

	if(ypn > .99e30)
	{
		qn=0.0;
      un=0.0;
      }
   else
   {
	   qn=0.50;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
      }
	
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

   for(int k=n-1;k>=1;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
	}



/************************************************************************************************************************************/
/*     rombint returns the integral from a to b of func using Romberg integration. The method converges provided that f(x) is 		*/
/*     continuous in (a,b). tol indicates the desired relative accuracy in the integral.															*/
/*		 The function is taken from : http://wise-obs.tau.ac.il/~barkana/InfallSub.c																	*/
/************************************************************************************************************************************/

double rombint(double (*func)(double), double a, double b, double tol)
{
	int MAXJ=5; 
  	int MAXITER=30,nint,i,j,k,jmax;
  	double g[MAXJ+1],h,gmax,g0,fourj,g1,error;
	
  	h=0.5*(b-a);
  	gmax=h*(func(a)+func(b));
  	g[0]=gmax;
  	nint=1;
  	error=1.e20;
  	for (i=0; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++)
  	{
		/*     Calculate next trapezoidal rule approximation to integral. */
      g0=0.0;
      for (k=1; k<=nint; k++)
		{
	 		 g0+=func(a+(k+k-1)*h);
	  /*	  printf("In rombint, %-15.8e\t%-15.8e\n",exp(a+(k+k-1)*h),g0);*/
			}
		   
      g0=0.5*g[0]+h*g0;
      h=0.5*h;
      nint*=2;
      jmax=(i<MAXJ)?i:MAXJ;
      fourj=1.;
      for (j=1; j<=jmax; j++)
	/*     Use Richardson extrapolation. */
		{
			fourj*=4.0;
		  	g1=g0+(g0-g[j-1])/(fourj-1.0);
		  	g[j-1]=g0;
		  	g0=g1;
			}
      /*      if (fabs(g0) > tol) */
//      printf("%d %e\n",i,fabs(g0));
	   if (fabs(g0) > 1.e-30) 
			error=1.0-gmax/g0;
      else 
			error=gmax; 
      gmax=g0;
      g[jmax]=g0;
    	}
    	
  	if ((i > MAXITER) && (fabs(error) > tol))
   	printf("rombint failed to converge; integral=%g, error=%g\n",g0,error);
  	return g0;
}

/*********************************************************************************************/
/*		Numerical recipe function.																					*/ 
/*    Used for sorting indexes																					*/
/*    http://www.koders.com/c/fid45CAF2ABADCA71599A411E7D9241EC88F4374C22.aspx?s=crc	      */
/*********************************************************************************************/

//#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define MSAN 7
#define NSTACK 50

void indexx(int n, double arr[], int indx[])
{
	int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	double a;

	istack=(int *)malloc(NSTACK*sizeof(int));

	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < MSAN) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=l;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1])
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) printf("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free(istack);
}
#undef MSAN
#undef NSTACK
#undef SWAP
//#undef NRANSI





/************************************************************************/
/*		This function is taken from the DVERK subroutine of CMBFAST	      */
/************************************************************************/


/************************************************************************/
/*                                                                      */
/*                                                                      */
/* if you discover any errors in this subroutine, please contact        */
/*                                                                      */
/*        Santanu Das                                           			*/
/*        Inter University Center for Astromomy and Astrophysics			*/
/*			 Post Bag 4																		*/
/*        Ganeshkhind																	*/
/*        Pune University Campus														*/
/*        Pune 411 007, India                                           */
/*                                                                      */
/*        phone: +91 20 2560 4516                                       */
/*                                                                      */
/*			 email:santanud@iucaa.ernet.in 											*/
/*					 sanjone@gmail.com													*/
/*                                                                      */
/* dverk is origionaly writtern by  KENNETH R JACKSON in the department */ 
/*	of computer science university of toronto. I was origionally written */
/*	in fortran 66. The code has been modified and written in C for using */
/*	in CMBANS by Santanu Das. The following note was also added by 		*/
/*	KENNETH R JACKSON 			                                    		*/
/*                                                                      */
/* the constants dwarf and rreb -- c(10) and c(11), respectively -- are */
/* set for a  vax  in  double  precision.  they  should  be  reset,  as */
/* described below, if this program is run on another machine.          */
/*                                                                      */
/* the c array is declared in this subroutine to have one element only, */
/* although  more  elements  are  referenced  in this subroutine.  this */
/* causes some compilers to issue warning messages.  there is,  though, */
/* no  error  provided  c is declared sufficiently large in the calling */
/* program, as described below.                                         */
/*                                                                      */
/* the following external statement  for  fcn  was  added  to  avoid  a */
/* warning  message  from  the  unix  f77 compiler.  the original dverk */
/* comments and code follow it.                                         */
/*                                                                      */
/************************************************************************/



/************************************************************************/
/*                                                                      */
/*     purpose - this is a runge-kutta  subroutine  based  on  verner's */
/* fifth and sixth order pair of formulas for finding approximations to */
/* the solution of  a  system  of  first  order  ordinary  differential */
/* equations  with  initial  conditions. it attempts to keep the global */
/* error proportional to  a  tolerance  specified  by  the  user.  (the */
/* proportionality  depends  on the kind of error control that is used, */
/* as well as the differential equation and the range of integration.)  */
/*                                                                      */
/*     various options are available to the user,  including  different */
/* kinds  of  error control, restrictions on step sizes, and interrupts */
/* which permit the user to examine the state of the  calculation  (and */
/* perhaps make modifications) during intermediate stages.              */
/*                                                                      */
/*     the program is efficient for non-stiff systems.  however, a good */
/* variable-order-adams  method  will probably be more efficient if the */
/* function evaluations are very costly.  such a method would  also  be */
/* more suitable if one wanted to obtain a large number of intermediate */
/* solution values by interpolation, as might be the case  for  example */
/* with graphical output.                                               */
/*                                                                      */
/*                                    hull-enright-jackson   1/10/76    */
/*                                                                      */
/************************************************************************/
/*                                                                      */
/*     use - the user must specify each of the following                */
/*                                                                      */
/*     n  number of equations                                           */
/*                                                                      */
/*   fcn  name of subroutine for evaluating functions - the  subroutine */
/*           itself must also be provided by the user - it should be of */
/*           the following form                                         */
/*              subroutine fcn(n, x, y, yprime)                         */
/*              integer n                                               */
/*              double precision x, y(n), yprime(n)                     */
/*                      *** etc ***                                     */
/*           and it should evaluate yprime, given n, x and y            */
/*                                                                      */
/*     x  independent variable - initial value supplied by user         */
/*                                                                      */
/*     y  dependent variable - initial values of components y(1), y(2), */
/*           ..., y(n) supplied by user                                 */
/*                                                                      */
/*  xend  value of x to which integration is to be carried out - it may */
/*           be less than the initial value of x                        */
/*                                                                      */
/*   tol  tolerance - the subroutine attempts to control a norm of  the */
/*           local  error  in  such  a  way  that  the  global error is */
/*           proportional to tol. in some problems there will be enough */
/*           damping  of  errors, as well as some cancellation, so that */
/*           the global error will be less than tol. alternatively, the */
/*           control   can   be  viewed  as  attempting  to  provide  a */
/*           calculated value of y at xend which is the exact  solution */
/*           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) */
/*           is proportional to tol.  (the norm  is  a  max  norm  with */
/*           weights  that  depend on the error control strategy chosen */
/*           by the user.  the default weight for the k-th component is */
/*           1/max(1,abs(y(k))),  which therefore provides a mixture of */
/*           absolute and relative error control.)                      */
/*                                                                      */
/*   ind  indicator - on initial entry ind must be set equal to  either */
/*           1  or  2. if the user does not wish to use any options, he */
/*           should set ind to 1 - all that remains for the user to  do */
/*           then  is  to  declare c and w, and to specify nw. the user */
/*           may also  select  various  options  on  initial  entry  by */
/*           setting ind = 2 and initializing the first 9 components of */
/*           c as described in the next section.  he may also  re-enter */
/*           the  subroutine  with ind = 3 as mentioned again below. in */
/*           any event, the subroutine returns with ind equal to        */
/*              3 after a normal return                                 */
/*              4, 5, or 6 after an interrupt (see options c(8), c(9))  */
/*              -1, -2, or -3 after an error condition (see below)      */
/*                                                                      */
/*       communications vector - the dimension must be greater than or  */
/*           equal to 24, unless option c(1) = 4 or 5 is used, in which */
/*           case the dimension must be greater than or equal to n+30   */
/*                                                                      */
/*    nw  first dimension of workspace w -  must  be  greater  than  or */
/*           equal to n                                                 */
/*                                                                      */
/*     w  workspace matrix - first dimension must be nw and second must */
/*           be greater than or equal to 9                              */
/*                                                                      */
/*     the subroutine  will  normally  return  with  ind  =  3,  having */
/* replaced the initial values of x and y with, respectively, the value */
/* of xend and an approximation to y at xend.  the  subroutine  can  be */
/* called  repeatedly  with new values of xend without having to change */
/* any other argument.  however, changes in tol, or any of the  options */
/* described below, may also be made on such a re-entry if desired.     */
/*                                                                      */
/*     three error returns are also possible, in which  case  x  and  y */
/* will be the most recently accepted values -                          */
/*     with ind = -3 the subroutine was unable  to  satisfy  the  error */
/*        requirement  with a particular step-size that is less than or */
/*        equal to hmin, which may mean that tol is too small           */
/*     with ind = -2 the value of hmin  is  greater  than  hmax,  which */
/*        probably  means  that the requested tol (which is used in the */
/*        calculation of hmin) is too small                             */
/*     with ind = -1 the allowed maximum number of fcn evaluations  has */
/*        been  exceeded,  but  this  can only occur if option c(7), as */
/*        described in the next section, has been used                  */
/*                                                                      */
/*     there are several circumstances that will cause the calculations */
/* to  be  terminated,  along with output of information that will help */
/* the user determine the cause of  the  trouble.  these  circumstances */
/* involve  entry with illegal or inconsistent values of the arguments, */
/* such as attempting a normal  re-entry  without  first  changing  the */
/* value of xend, or attempting to re-enter with ind less than zero.    */
/*                                                                      */
/************************************************************************/
/*                                                                      */
/*     options - if the subroutine is entered with ind = 1, the first 9 */
/* components of the communications vector are initialized to zero, and */
/* the subroutine uses only default values  for  each  option.  if  the */
/* subroutine  is  entered  with ind = 2, the user must specify each of */
/* these 9 components - normally he would first set them all  to  zero, */
/* and  then  make  non-zero  those  that  correspond to the particular */
/* options he wishes to select. in any event, options may be changed on */
/* re-entry  to  the  subroutine  -  but if the user changes any of the */
/* options, or tol, in the course of a calculation he should be careful */
/* about  how  such changes affect the subroutine - it may be better to */
/* restart with ind = 1 or 2. (components 10 to 24 of c are used by the */
/* program  -  the information is available to the user, but should not */
/* normally be changed by him.)                                         */
/*                                                                      */
/*  c(1)  error control indicator - the norm of the local error is  the */
/*           max  norm  of  the  weighted  error  estimate  vector, the */
/*           weights being determined according to the value of c(1) -  */
/*              if c(1)=1 the weights are 1 (absolute error control)    */
/*              if c(1)=2 the weights are 1/abs(y(k))  (relative  error */
/*                 control)                                             */
/*              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) */
/*                 (relative  error  control,  unless abs(y(k)) is less */
/*                 than the floor value, abs(c(2)) )                    */
/*              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) */
/*                 (here individual floor values are used)              */
/*              if c(1)=5 the weights are 1/abs(c(k+30))                */
/*              for all other values of c(1), including  c(1) = 0,  the */
/*                 default  values  of  the  weights  are  taken  to be */
/*                 1/max(1,abs(y(k))), as mentioned earlier             */
/*           (in the two cases c(1) = 4 or 5 the user must declare  the */
/*           dimension of c to be at least n+30 and must initialize the */
/*           components c(31), c(32), ..., c(n+30).)                    */
/*                                                                      */
/*  c(2)  floor value - used when the indicator c(1) has the value 3    */
/*                                                                      */
/*  c(3)  hmin specification - if not zero, the subroutine chooses hmin */
/*           to be abs(c(3)) - otherwise it uses the default value      */
/*              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     */
/*           where dwarf is a very small positive  machine  number  and */
/*           rreb is the relative roundoff error bound                  */
/*                                                                      */
/*  c(4)  hstart specification - if not zero, the subroutine  will  use */
/*           an  initial  hmag equal to abs(c(4)), except of course for */
/*           the restrictions imposed by hmin and hmax  -  otherwise it */
/*          uses the default value of hmax*(tol)**(1/6)                 */
/*                                                                      */
/*  c(5)  scale specification - this is intended to be a measure of the */
/*           scale of the problem - larger values of scale tend to make */
/*           the method more reliable, first  by  possibly  restricting */
/*           hmax  (as  described  below) and second, by tightening the */
/*           acceptance requirement - if c(5) is zero, a default  value */
/*           of  1  is  used.  for  linear  homogeneous  problems  with */
/*          constant coefficients, an appropriate value for scale is a  */
/*           norm  of  the  associated  matrix.  for other problems, an */
/*           approximation to  an  average  value  of  a  norm  of  the */
/*           jacobian along the trajectory may be appropriate           */
/*                                                                      */
/*  c(6)  hmax specification - four cases are possible                  */
/*          if c(6).ne.0 and c(5).ne.0, hmax is taken to be             */
/*              min(abs(c(6)),2/abs(c(5)))                              */
/*           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) */
/*           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            */
/*              2/abs(c(5))                                             */
/*           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value */
/*              of 2                                                    */
/*                                                                      */
/*  c(7)  maximum number of function evaluations  -  if  not  zero,  an */
/*           error  return with ind = -1 will be caused when the number */
/*           of function evaluations exceeds abs(c(7))                  */
/*                                                                      */
/*  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will */
/*           interrupt   the  calculations  after  it  has  chosen  its */
/*           preliminary value of hmag, and just before choosing htrial */
/*           and  xtrial  in  preparation for taking a step (htrial may */
/*           differ from hmag in sign, and may  require  adjustment  if */
/*           xend  is  near) - the subroutine returns with ind = 4, and */
/*           will resume calculation at the point  of  interruption  if */
/*           re-entered with ind = 4                                    */
/*                                                                      */
/*  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will */
/*           interrupt   the  calculations  immediately  after  it  has */
/*           decided whether or not to accept the result  of  the  most */
/*           recent  trial step, with ind = 5 if it plans to accept, or */
/*           ind = 6 if it plans to reject -  y(*)  is  the  previously */
/*           accepted  result, while w(*,9) is the newly computed trial */
/*           value, and w(*,2) is the unweighted error estimate vector. */
/*           the  subroutine  will  resume calculations at the point of */
/*           interruption on re-entry with ind = 5 or 6. (the user  may */
/*           change ind in this case if he wishes, for example to force */
/*           acceptance of a step that would otherwise be rejected,  or */
/*           vice versa. he can also restart with ind = 1 or 2.)        */
/*                                                                      */
/************************************************************************/
/*                                                                      */
/*  summary of the components of the communications vector              */
/*                                                                      */
/*     prescribed at the option       determined by the program         */
/*           of the user                                                */
/*                                                                      */
/*                                    c(10) rreb(rel roundoff err bnd)  */
/*     c(1) error control indicator   c(11) dwarf (very small mach no)  */
/*     c(2) floor value               c(12) weighted norm y             */
/*     c(3) hmin specification        c(13) hmin                        */
/*     c(4) hstart specification      c(14) hmag                        */
/*     c(5) scale specification       c(15) scale                       */
/*     c(6) hmax specification        c(16) hmax                        */
/*     c(7) max no of fcn evals       c(17) xtrial                      */
/*     c(8) interrupt no 1            c(18) htrial                      */
/*     c(9) interrupt no 2            c(19) est                         */
/*                                    c(20) previous xend               */
/*                                    c(21) flag for xend               */
/*                                    c(22) no of successful steps      */
/*                                    c(23) no of successive failures   */
/*                                    c(24) no of fcn evals             */
/*                                                                      */
/*  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        */
/*                                                                      */
/************************************************************************/
/*                                                                      */
/*  an overview of the program                                          */
/*                                                                      */
/*     begin initialization, parameter checking, interrupt re-entries   */
/*  ......abort if ind out of range 1 to 6                              */
/* .     cases - initial entry, normal re-entry, interrupt re-entries   */
/*  .     case 1 - initial entry (ind .eq. 1 or 2)                      */
/*  v........abort if n.gt.nw or tol.le.0                               */
/*  .        if initial entry without options (ind .eq. 1)              */
/*  .           set c(1) to c(9) equal to zero                          */
/*  .        else initial entry with options (ind .eq. 2)               */
/*  .           make c(1) to c(9) non-negative                          */
/*  .           make floor values non-negative if they are to be used   */
/*  .        end if                                                     */
/*  .        initialize rreb, dwarf, prev xend, flag, counts            */
/*  .     case 2 - normal re-entry (ind .eq. 3)                         */
/*  .........abort if xend reached, and either x changed or xend not    */
/*  .        re-initialize flag                                         */
/*  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    */
/*  v        transfer control to the appropriate re-entry point.......  */
/*  .     end cases                                                  .  */
/*  .  end initialization, etc.                                      .  */
/*  .                                                                v  */
/*  .  loop through the following 4 stages, once for each trial step .  */
/*  .     stage 1 - prepare                                          .  */
/************error return (with ind=-1) if no of fcn evals too great .  */
/*  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  */
/*  .        calc hmin, scale, hmax                                  .  */
/************error return (with ind=-2) if hmin .gt. hmax            .  */
/*  .        calc preliminary hmag                                   .  */
/************interrupt no 1 (with ind=4) if requested.......re-entry.v  */
/*  .        calc hmag, xtrial and htrial                            .  */
/*  .     end stage 1                                                .  */
/*  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  */
/*  .     stage 3 - calc the error estimate                          .  */
/*  .     stage 4 - make decisions                                   .  */
/*  .        set ind=5 if step acceptable, else set ind=6            .  */
/************interrupt no 2 if requested....................re-entry.v  */
/*  .        if step accepted (ind .eq. 5)                              */
/*  .           update x, y from xtrial, ytrial                         */
/*  .           add 1 to no of successful steps                         */
/*  .           set no of successive failures to zero                   */
/***************return(with ind=3, xend saved, flag set) if x .eq. xend */
/*  .        else step not accepted (ind .eq. 6)                        */
/*  .           add 1 to no of successive failures                      */
/***************error return (with ind=-3) if hmag .le. hmin            */
/*  .        end if                                                     */
/*  .     end stage 4                                                   */
/*  .  end loop                                                         */
/*  .                                                                   */
/*  begin abort action                                                  */
/*     output appropriate  message  about  stopping  the  calculations, */
/*       along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend,  */
/*       previous xend,  no of  successful  steps,  no  of  successive  */
/*       failures, no of fcn evals, and the components of y             */
/*     stop                                                             */
/*  end abort action                                                    */
/*                                                                      */
/************************************************************************/


void testflg(int nnn, int mm,double cc[],double ww[][nvar0+1])
{
	printf("c:\n");
	for(int i=1;i<=24;i+=4)
		printf("[%d] = %e,[%d] = %e,[%d] = %e,[%d] = %e\n",i,cc[i],i+1,cc[i+1],i+2,cc[i+2],i+3,cc[i+3]);
	printf("W:\n");	
	for(int i=1;i<=nnn;i+=4)
		printf("[%d] = %e,[%d] = %e,[%d] = %e,[%d] = %e\n",i,ww[mm][i],i+1,ww[mm][i+1],i+2,ww[mm][i+2],i+3,ww[mm][i+3]);
	}
	
void errmsg(int n,double  x,double y[],double xend,double tol,int  ind,double c[],int nw)
{	
	printf("Some error occure,\n");
	printf("computation stopped in dverk with the following values :-\n");
	printf("ind = %d\t tol  =%e\t x =%e\t n =%d\n",ind,tol,x,n);
	printf("hmin = %e\t xend = %e\tnw  =%d\t hmax = %e\n",c[13],xend,nw,c[16]);
	printf("prev xend =%e",c[20]);
   printf("no of successful steps    =%e\n",c[22]);
   printf("no of successive failures =%e\n",c[23]);
   printf("no of function evals      =%e\n",c[24]);
   
   printf("the components of y are");
   for(int i=1;i<=n;i++)
   	printf("y[%d] = %e",i,y[i]);
	
	exit(1);
	}
	
	
double dverk (int nsann,int n,void (*fcn)(int, double, double[],double[]),double x,double y[],double xend,double tol,int  ind,double c[],int nw,double w[][nvar0+1])
{
//	int k;
   double temp;
   int run=1;
	int santest=0;
	
	if(xend>=206.1075)
		santest=1;
	/******************************************************************************/
	/*      begin initialization, parameter checking, interrupt re-entries 			*/
	/******************************************************************************/
	
   if (ind<1 || ind>6)																//  ......abort if ind out of range 1 to 6
   	errmsg(n,x,y,xend,tol,ind,c,nw);
   
   																						// cases - initial entry, normal re-entry, interrupt re-entries
   if(ind == 4) goto A111;
   else if(ind ==5 || ind == 6) goto A222;
   
	if(ind==1 || ind ==2)															//  .........case 1 - initial entry (ind .eq. 1 or 2)
	{
 		if(n>nw || tol<=0.0)															//  .........abort if n.gt.nw or tol.le.0
      	errmsg(n,x,y,xend,tol,ind,c,nw);
      	
     	if (ind==1)
      {
	      for(int k=1; k<=9; k++)													//  initial entry without options (ind .eq. 1)
      		c[k] = 0.0;		
//      	printf("Hi\n");															//  set c(1) to c(9) equal to 0
			}
				
		else
		{
  			for(int k=1; k<=9; k++)													//  initial entry with options (ind .eq. 2)
      		c[k] = fabs(c[k]);														//  make c(1) to c(9) non-negative
   		
      	if (c[1]==4.0 || c[1]==5.0)											//	 make floor values non-negative if they are to be used
     		{
	     		for(int k=1;k<=n;k++)
         	{
	         	c[k+30] = fabs(c[k+30]);
   	      	}
				}
			}
			
//   	printf("Pass 2\n");		
  		c[10] = pow(2.0,-56);//1.38777878e-17;   													// initialize rreb, dwarf, prev xend, flag, counts
      c[11] = 1.0e-35;																
//    printf("%e %e",c[10],c[11]);
      c[20] = x;																		//  set previous xend initially to initial value of x
      
      for(int k=21; k<=24; k++)
		{
			c[k] = 0.0;
   		}
		}
																							//  case 2 - normal re-entry (ind .eq. 3)
																							// .........abort if xend reached, and either x changed or xend not
	if(ind ==3)
	{		
		if (c[21] != 0.0 && (x!=c[20] || xend==c[20]))
			errmsg(n,x,y,xend,tol,ind,c,nw);																	//  re-initialize flag
      c[21] = 0.0;
      }

	//  case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
	//  transfer control to the appropriate re-entry point..........
	//  this has already been handled by the computed go to. end cases
	
	//  end initialization, etc.


	/*********************************************************************************/
	/*       loop through the following 4 stages, once for each trial  step  			*/
	/*       until the occurrence of one of the following                    			*/	
	/*          (a) the normal return (with ind .eq. 3) on reaching xend in  			*/
	/*              stage 4                                                  			*/
	/*          (b) an error return (with ind .lt. 0) in stage 1 or stage 4  			*/
	/*          (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if  			*/
	/*              requested, in stage 1 or stage 4                         			*/
	/*********************************************************************************/

	while(1)
	{
//	printf("%d x = %e xend = %e c[14] = %e c[1] = %e c[13] = %e c[7] = %e \nc[24] = %e c[3] = %e c[12] = %e c[15] = %e c[16] = %e c[6] = %e c[23] = %e  c[5] = %e",run,x,xend,c[14],c[1],c[13],c[7],c[24],c[3],c[12],c[15],c[16],c[6],c[23],c[5]);
//	printf("\nind = %d c[14] = %e",ind,c[14]);
		/*********************************************************************************/
		/*          stage 1 - prepare - do calculations of  hmin,  hmax,  etc.,  			*/	
		/*          and some parameter  checking,  and  end  up  with  suitable  			*/
		/*          values of hmag, xtrial and htrial in preparation for taking  			*/	
		/*          an integration step.                                         			*/
		/*********************************************************************************/   

		if (c[7]!=0.0 && c[24]>=c[7])												//***********error return (with ind=-1) if no of fcn evals too great
		{	
			ind = -1;
   	   return x;
   	   }
      		
   	if (ind != 6)																	//  calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
   	{
		   fcn(n, x, y, w[1]);
   	   c[24] = c[24] + 1.0;
   	   }
   	

				   
		c[13] = c[3];																	//  calculate hmin - use default unless value prescribed
//		printf("1 C[13] = %e C[16] = %e\n",c[13],c[16]);
			
		if (c[3] == 0.0)
		{
			//  calculate default value of hmin
			//  first calculate weighted norm y - c(12) - as specified
			//  by the error control indicator c(1)
   	   	
   	   temp = 0.0;
//   	   printf("C[1] = %e\t",c[1]);
		 	if (c[1] == 1.0)
  		   {																	  		   //	 absolute error control - weights are 1
      		for(int k=1; k<=n; k++)
      			temp = max(temp,fabs(y[k]));
	      	c[12] = temp;
				}
					
			else if (c[1] == 2.0)
			{																				//  relative error control - weights are 1/abs(y(k)) so
	         c[12] = 1.0;															//  weighted norm y is 1
	         }
			
			else if (c[1] == 3.0)
			{																				//	weights are 1/max(c(2),abs(y(k)))
	         for(int k=1; k<=n; k++)
		         temp = max(temp, fabs(y[k])/c[2]);
				c[12] = min(temp, 1.0);
	         }
  	       
			else if (c[1] == 4.0)
			{																				//	weights are 1/max(c(k+30),abs(y(k)))
  	   	 	for(int k=1; k<=n; k++)
	   	   	temp = max(temp, fabs(y[k])/c[k+30]);
      		c[12] = min(temp, 1.0);
      		}

			else if (c[1] == 5.0)
			{																				// weights are 1/c(k+30)
   	     	for(int k=1; k<=n; k++)
   	      	temp = max(temp, fabs(y[k]/c[k+30]));
				c[12] = temp;
   	      }
					
			else																			// default case - weights are 1/max(1,abs(y(k)))
			{
				for(int k=1;k<=n;k++)
   	   		temp = max(temp, fabs(y[k]));
			 	c[12] = min(temp, 1.0);
 				} 		
// 			printf("%e %e %e %e %e\n",c[11],c[10],c[12],x,tol);	
 			
			c[13] = 10.0*max(c[11],c[10]*max(c[12]/tol,fabs(x)));
	  		}

//		printf("\nc[14] new %e",c[14]);
	   c[15] = c[5];																	// calculate scale - use default unless value prescribed
//		printf("2 C[13] = %e C[16] = %e\n",c[13],c[16]);
			   		
	   if(c[5] == 0.0) 
			c[15] = 1.0;
																									// calculate hmax - consider 4 cases
		if(c[6] != 0.0 && c[5]!=0.0)												// case 1 both hmax and scale prescribed
	    	c[16] = min(c[6], 2.0/c[5]);
						      
		if (c[6] != 0.0 && c[5]==0.0) 											// case 2 - hmax prescribed, but scale not
			c[16] = c[6];
						      
		if (c[6] == 0.0 && c[5] != 0.0)									 		// case 3 - hmax not prescribed, but scale is
			c[16] = 2.0/c[5];
	


		if (c[6] == 0.0 && c[5] == 0.0) 											// case 4 - neither hmax nor scale is provided
			c[16] = 2.0;
			


	 	if (c[13] > c[16])															//***********error return (with ind=-2) if hmin .gt. hmax
	   {
		   ind = -2;
//		   printf("Pass - this");
	      return x;
	  		}
//	printf("\n%e %e %d\n",c[14],c[4],ind);
//		printf("ind = %d",ind);				
	   if (ind <= 2)																	// calculate preliminary hmag - consider 3 cases
	   {																					// case 1 - initial entry - use prescribed value of hstart, if
	      c[14] = c[4];																// any, else default
	      if (c[4] == 0.0) 
	      	c[14] = c[16]*pow(tol,(1.0/6.0));
//	      goto A444;	
	 	   }
		else if(c[23]<=1.0)
	  	{																					// case 2 - after a successful step, or at most  one  failure,
	      temp = 2.0*c[14];															// use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
	      if (tol < pow((2.0/.90),6)*c[19])		
	      {							// overflow. then avoid reduction by more than half.
	      	temp = .90*pow((tol/c[19]),(1.0/6.0))*c[14];
	     		}
	      c[14] = max(temp,0.5*c[14]);
//     		goto A444;
			}
	   		

//		else if(tol >= pow((2.0/.90),6)*c[19])	/// Checking important
//printf("\n%e\n",c[14]);
		else 
      	c[14] = 0.50*c[14];														// case 3 - after two or more successive failures
//printf("\n%e\n",c[14]);
//A444:	
	   c[14] = min(c[14],c[16]);													// check against hmax
	
	   c[14] = max(c[14], c[13]);													// check against hmin
//		printf("%e %e %e %e %e\n",c[6],c[5],c[16],c[14],c[13]);
//		exit(1);		
		
		
	   if (c[8] != 0.0)																//***********interrupt no 1 (with ind=4) if requested
	   {
		   ind = 4;
//			printf("A23");	   		   
	      return x;																		// resume here on re-entry with ind .eq. 4   ........re-entry..
			}
//		printf("Hellow -3");	   

//		printf("c14 = %e diff = %e\n",c[14],fabs(xend - x));
A111: 
//		printf("\nc[14] new %e %d",c[14],ind);
	   if (c[14] < fabs(xend - x))												// calculate hmag, xtrial - depending on preliminary hmag, xend
		{	
//			printf(" Error here 1  %e ",c[17]);
//			printf("Hi %e %e\n",c[14],.50*fabs(xend - x));									
	      c[14] = min(c[14], .50*fabs(xend - x));								// do not step more than half way to xend
	      c[17] = x + (((xend - x)>0.0)?fabs(c[14]):-fabs(c[14]));
	      }
	      
		else																				// hit xend exactly
		{		
//			printf(" Error here 2  %e ",c[17]);		
//			printf("Hi1");
			c[14] = fabs(xend-x);	
	      c[17] = xend;
  		   }

	   c[18] = c[17] - x;															// calculate htrial

 //  		if(run==2)
//	  		{
//		  		varify(n,x,y,xend,tol,ind,c,nw);
//				exit(1);
//				}														
		// end stage 1
	
		/*        ***************************************************************			*/
		/*        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *			*/
		/*        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *			*/
		/*        * stage 3. w(*,9) is temporary storage until finally it holds *			*/
		/*        * ytrial.                                                     *			*/
		/*        ***************************************************************			*/
	
	   temp = c[18]/1398169080000.0;
	   
//		printf("\ntemp = %e x = %e c18/6 = %e\n",temp,x,c[18]/6);

		
		for(int k=1; k<=n; k++)
	   	w[9][k] = y[k] + temp*w[1][k]*233028180000.0;

//		testflg(n,9,c,w);	
		
 	   fcn(n, x + c[18]/6.0, w[9], w[2]);
//		printf("\nnext time = %e",x + c[18]/6.0);
//		testflg(n,2,c,w);		   
//		testflg(n,9,c,w);		   
//		exit(1);
//		for(int ih=1;ih<=9;ih+=3)
//		{
//			printf("[%d]= %.8e [%d]= %.8e [%d]= %.8e\n",ih,w[ih][1],ih+1,w[ih+1][1],ih+2,w[ih+2][1]);
//			}	

	   for(int k=1; k<=n; k++)
	   	w[9][k] = y[k] + temp*(w[1][k]*74569017600.0 + w[2][k]*298276070400.0);
	  
	   fcn(n, x + c[18]*(4.0/15.0), w[9], w[3]);

		 	
	   for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(w[1][k]*1165140900000.0 - w[2][k]*3728450880000.0 + w[3][k]*3495422700000.0);	

	   fcn(n, x + c[18]*(2.0/3.0), w[9], w[4]);
	
	   for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(- w[1][k]*3604654659375.0 + w[2][k]*12816549900000.0 - w[3][k]*9284716546875.0 + w[4][k]*1237962206250.0);
	               
	   fcn(n, x + c[18]*(5.0/6.0), w[9], w[5]);
	
	   for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(w[1][k]*3355605792000.0 - w[2][k]*11185352640000.0 + w[3][k]*9172628850000.0 - w[4][k]*427218330000.0 + w[5][k]*482505408000.0);
	
	   fcn(n, x + c[18], w[9], w[6]);

		
	   for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(-w[1][k]*770204740536.0 + w[2][k]*2311639545600.0 - w[3][k]*1322092233000.0 - w[4][k]*453006781920.0 + w[5][k]*326875481856.0 );
	
	   fcn(n, x + c[18]/15.0, w[9], w[7]);

		
	   for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(w[1][k]*2845924389000.0 - w[2][k]*9754668000000.0 + w[3][k]*7897110375000.0 - w[4][k]*192082660000.0  + w[5][k]*400298976000.0 + w[7][k]*201586000000.0 );
	  
		fcn(n, x + c[18], w[9], w[8]);	

		//  calculate ytrial, the extrapolated approximation and store
		//  in w(*,9)
		for(int k=1;k<=n;k++)
	   	w[9][k] = y[k] + temp*(w[1][k]*104862681000.0 + w[3][k]*545186250000.0 + w[4][k]*446637345000.0 + w[5][k]*188806464000.0 + w[7][k]*15076875000.0 + w[8][k]*97599465000.0 );
	
		//  add 7 to the no of fcn evals
//		if(nsann==12)
//		{
//		for(int sanint=1;sanint<=41;sanint++)
//		{
//				for(int intsan=1;intsan<=9;intsan++)
//					printf("%.8e    ",w[intsan][sanint]);
//				printf("\n");
//				}	 
//			exit(1);
//			}
//		printf("\n\n\n\n");				 		
		
//		testflg(n,1,c,w);		   

	   c[24] = c[24] + 7.0;

		//        end stage 2

	
		/*        ***************************************************************			*/
		/*        * stage 3 - calculate the error estimate est. first calculate *			*/
		/*        * the  unweighted  absolute  error  estimate vector (per unit *			*/	
		/*        * step) for the unextrapolated approximation and store it  in *			*/
		/*        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *			*/
		/*        * specified by the error  control  indicator  c(1).  finally, *			*/
		/*        * modify  this result to produce est, the error estimate (per *			*/
		/*        * unit step) for the extrapolated approximation ytrial.       *			*/
		/*        ***************************************************************			*/
	
		//           calculate the unweighted absolute error estimate vector
	   for(int k=1;k<=n;k++)
	   	w[2][k] = (w[1][k]*(8738556750.0/1398169080000.0) + w[3][k]*(9735468750.0/1398169080000.0) - w[4][k]*(9709507500.0/1398169080000.0) + w[5][k]*(8582112000.0/1398169080000.0) + w[6][k]*(95329710000.0/1398169080000.0) - w[7][k]*(15076875000.0/1398169080000.0) - w[8][k]*(97599465000.0/1398169080000.0));

		//           calculate the weighted max norm of w(*,2) as specified by
		//           the error control indicator c(1)
	   temp = 0.0;
//		printf("Hellow -2");	   

	   if(c[1] == 1.0)
	   {																					// absolute error control
	      for(int k=1;k<= n;k++)
	      	temp = max(temp,fabs(w[2][k]));
	  		 }
	  	
		else if (c[1] == 2.0)
	  	{																					// relative error control
	      for(int k=1;k<=n;k++)
	      	temp = max(temp, fabs(w[2][k]/y[k]));
			}
	   
		else if (c[1] == 3.0) 
	  	{																					// weights are 1/max(c(2),abs(y(k)))
	      for(int k=1;k<=n;k++)
	      	temp = max(temp, fabs(w[2][k])/ max(c[2], fabs(y[k])));
	      }
	   
	   else if(c[1] == 4.0)
	  	{																					// weights are 1/max(c(k+30),abs(y(k)))
	      for(int k=1;k<=n;k++)
	      	temp = max(temp, fabs(w[k][2])/max(c[k+30],fabs(y[k])));
			}
	
		else if(c[1] == 5.0)
	  	{																					// weights are 1/c(k+30)
	      for(int k=1;k<=n;k++)
	      	temp = max(temp, fabs(w[2][k]/c[k+30]));
			}
	
		else
		{		
																	   	// default case - weights are 1/max(1,abs(y(k)))
	      for(int k=1;k<=n;k++)
	      {
//				printf("\nHi,k = %d temp = %e w[2][k] = %e y[k] = %e\n",k,temp,w[2][k],y[k]);	      
	      	temp = max(temp, fabs(w[2][k])/max(1.0, fabs(y[k])));
	      	}
	      }

	   //
		//           calculate est - (the weighted max norm of w(*,2))*hmag*scale
		//              - est is intended to be a measure of the error  per  unit
		//              step in ytrial
	   c[19] = temp*c[14]*c[15];

	
//		printf("\nc[19] = %e temp = %e c14 = %e c15 = %e\n",c[19],temp,c[14],c[15]);
		
		//        end stage 3
		
		/*        ***************************************************************			*/
		/*        * stage 4 - make decisions.                                   *			*/
		/*        ***************************************************************			*/
	
		//           set ind=5 if step acceptable, else set ind=6
	   
	   ind = 5;
	   
	   if (c[19]>tol) 
	   	ind = 6;
	
		//***********interrupt no 2 if requested
		if(c[9] != 0.0) 
		{
//			printf("Pass A1");
	   	return x;
			}
		// resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
		//
A222:		
//		printf("Hellow -1");
//		printf("%d\n",ind);
		if(ind != 6) 
		{
			
			// step accepted (ind .eq. 5), so update x, y from xtrial,
			// ytrial, add 1 to the no of successful steps, and set
			// the no of successive failures to zero
	      x = c[17];
			for(int k=1;k<=n;k++)
	      	y[k] = w[9][k];
	
			c[22] = c[22] + 1.0;
			c[23] = 0.0;
			//**************return(with ind=3, xend saved, flag set) if x .eq. xend
	      if (x == xend)
	      {
		      ind = 3;
	         c[20] = xend;
	         c[21] = 1.0;
//				printf("Pass A2");
	         return x;
	         }
			}
	   else
	   {
			// step not accepted (ind .eq. 6), so add 1 to the no of
			// successive failures
	      c[23] = c[23] + 1.0;
	      //**************error return (with ind=-3) if hmag .le. hmin
	      if(c[14] <= c[13])
	      {
		      ind = -3;
//				printf("Pass A3");
	         return x;
				}
	  		}

			
		// end stage 4
//		printf("Hellow");
//		if(santest==1)
//		{	
//			printf("\nrun = %d tauend = %e y = \n\n",run,xend);	
/*			for(int abcd=1;abcd<=n;abcd++)
			   printf("%e\t",y[abcd]);
			printf("\n Now W \n");
			for(int bcd=1;bcd<=9;bcd++)
			{
			printf("\nW%d",bcd);	
			for(int abcd=1;abcd<=n;abcd++)
			   printf("%e\t",w[bcd][abcd]);
			 }     
			printf("\n");*/
//		if(run==100)exit(1); 
//			}
		run++;	  

		}	//     end loop
	
	//  begin abort action
	}


#endif