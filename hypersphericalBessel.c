/* 
 * Calculation of hyperspherical Bessel functions, including
 * spherical Bessel functions in the flat-space case.
 *
 * The routine phi_recursive uses recursion relations to calculate
 * the functions, which is accurate but relatively slow.
 * CHECK ACCURACY FOR BETA < 1
 *
 * The routine phi_WKB uses Langer's formula for a
 * uniform first-order asymptotic approximation in the open and
 * flat cases. This approximation is EXCELLENT for all l >= 3.
 * 
 * The routine qintegral calculates the closed-form answer
 * to the eikonal integral used in the WKB approximation.
 *  
 * The routine qintegral_simp calculates the eikonal integral used in the
 * WKB approximation via a straightforward Simpson's Rule integration.
 * ** For test purposes only **
 *
 * The routine airy_ai returns the Airy function Ai(x) of the argument 
 * passed. It employs a Pade-type approximation away from zero and 
 * a Taylor expansion around zero. Highly accurate.
 *
 * The routines polevl and p1evl are auxiliary polynomial
 * evaluation routines used in the airy function calculation.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>




// AIRY FUNCTION 


/*	polevl.c
 *	p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 *  The function p1evl() assumes that coef[N] = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


double polevl(double x, double coef[], int N )
{
double ans;
int i;
double *p;

p = coef;
ans = *p++;
i = N;

do
	ans = ans * x  +  *p++;
while( --i );

return( ans );
}

/*					      
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl(double x, double coef[], int N )
{
double ans;
double *p;
int i;

p = coef;
ans = x + *p++;
i = N-1;

do
	ans = ans * x  + *p++;
while( --i );

return( ans );
}



/***********************************************************************
 *                                                                     *
 *	Airy function                                                  *
 *                                                                     *
 * Modified from original routine by Stephen Moshier, available        *
 * as part of the Cephes library at www.netlib.com                     *
 * Modifications: eliminates calculation of Bi(x), Ai'(x), Bi'(x)      *
 *                                                                     *
 ***********************************************************************/

/*
 * SYNOPSIS:
 *
 * double x, ai, aip, bi, bip;
 * int airy();
 *
 * airy( x, _&ai, _&aip, _&bi, _&bip );
 *
 *
 *
 * DESCRIPTION:
 *
 * Solution of the differential equation
 *
 *	y"(x) = xy.
 *
 * The function returns the two independent solutions Ai, Bi
 * and their first derivatives Ai'(x), Bi'(x).
 *
 * Evaluation is by power series summation for small x,
 * by rational minimax approximations for large x.
 *
 *
 *
 * ACCURACY:
 * Error criterion is absolute when function <= 1, relative
 * when function > 1, except * denotes relative error criterion.
 * For large negative x, the absolute error increases as x^1.5.
 * For large positive x, the relative error increases as x^1.5.
 *
 * Arithmetic  domain   function  # trials      peak         rms
 * IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16
 * IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*
 * IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16
 * IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*
 * IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16
 * IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16
 * DEC         -10, 0     Ai         5000       1.7e-16     2.8e-17
 * DEC           0, 10    Ai         5000       2.1e-15*    1.7e-16*
 * DEC         -10, 0     Ai'        5000       4.7e-16     7.8e-17
 * DEC           0, 10    Ai'       12000       1.8e-15*    1.5e-16*
 * DEC         -10, 10    Bi        10000       5.5e-16     6.8e-17
 * DEC         -10, 10    Bi'        7000       5.3e-16     8.7e-17
 *
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
/* Download from: http://www.moshier.net/#Cephes  */

/**************************************************/
/*  Calculation of the Airy function              */
/**************************************************/

static double c1 = 0.35502805388781723926;
static double c2 = 0.258819403792806798405;
static double sqrt3 = 1.732050807568877293527;
static double sqpii = 5.64189583547756286948E-1;

#define MAXAIRY 25.77
#define ACC 1.e-8
#define PI 3.1415926536

static double AN[8] = {
  3.46538101525629032477E-1,
  1.20075952739645805542E1,
  7.62796053615234516538E1,
  1.68089224934630576269E2,
  1.59756391350164413639E2,
  7.05360906840444183113E1,
  1.40264691163389668864E1,
  9.99999999999999995305E-1,
};

static double AD[8] = {
  5.67594532638770212846E-1,
  1.47562562584847203173E1,
  8.45138970141474626562E1,
  1.77318088145400459522E2,
  1.64234692871529701831E2,
  7.14778400825575695274E1,
  1.40959135607834029598E1,
  1.00000000000000000470E0,
};

static double AFN[9] = {
-1.31696323418331795333E-1,
-6.26456544431912369773E-1,
-6.93158036036933542233E-1,
-2.79779981545119124951E-1,
-4.91900132609500318020E-2,
-4.06265923594885404393E-3,
-1.59276496239262096340E-4,
-2.77649108155232920844E-6,
-1.67787698489114633780E-8,
};

static double AFD[9] = {
/* 1.00000000000000000000E0,*/
 1.33560420706553243746E1,
 3.26825032795224613948E1,
 2.67367040941499554804E1,
 9.18707402907259625840E0,
 1.47529146771666414581E0,
 1.15687173795188044134E-1,
 4.40291641615211203805E-3,
 7.54720348287414296618E-5,
 4.51850092970580378464E-7,
};

static double AGN[11] = {
  1.97339932091685679179E-2,
  3.91103029615688277255E-1,
  1.06579897599595591108E0,
  9.39169229816650230044E-1,
  3.51465656105547619242E-1,
  6.33888919628925490927E-2,
  5.85804113048388458567E-3,
  2.82851600836737019778E-4,
  6.98793669997260967291E-6,
  8.11789239554389293311E-8,
  3.41551784765923618484E-10,
};

static double AGD[10] = {
/*  1.00000000000000000000E0,*/
  9.30892908077441974853E0,
  1.98352928718312140417E1,
  1.55646628932864612953E1,
  5.47686069422975497931E0,
  9.54293611618961883998E-1,
  8.64580826352392193095E-2,
  4.12656523824222607191E-3,
  1.01259085116509135510E-4,
  1.17166733214413521882E-6,
  4.91834570062930015649E-9,
};


double airy_ai(double x)
{
  double z, zz, t, f, g, uf, ug, k, zeta, theta;
  double ai;

  /* Exponentially tiny for large enough argument
   */
  if( x > MAXAIRY ) return 0.;

  /* Pade fit for large negative arguments
   */
  if( x < -2.09 ){
    t = sqrt(-x);
    zeta = -2.0 * x * t / 3.0;
    t = sqrt(t);
    k = sqpii / t;
    z = 1.0/zeta;
    zz = z * z;
    uf = 1.0 + zz * polevl( zz, AFN, 8 ) / p1evl( zz, AFD, 9 );
    ug = z * polevl( zz, AGN, 10 ) / p1evl( zz, AGD, 10 );
    theta = zeta + 0.25 * PI;
    f = sin( theta );
    g = cos( theta );
    ai = k * (f * uf - g * ug);

    return(ai);
  }

  /* Pade fit for large positive arguments
   */
  if( x >= 2.09 ){	/* cbrt(9) */
    t = sqrt(x);
    zeta = 2.0 * x * t / 3.0;
    g = exp( zeta );
    t = sqrt(t);
    k = 2.0 * t * g;
    z = 1.0/zeta;
    f = polevl( z, AN, 7 ) / polevl( z, AD, 7 );
    ai = sqpii * f / k;
    return(ai);
  }


  /* Taylor series for region around x=0
   */

  f = 1.0;
  g = x;
  t = 1.0;
  uf = 1.0;
  ug = x;
  k = 1.0;
  z = x * x * x;
  while( t > ACC ){
    uf *= z;
    k += 1.0;
    uf /=k;
    ug *= z;
    k += 1.0;
    ug /=k;
    uf /=k;
    f += uf;
    k += 1.0;
    ug /=k;
    g += ug;
    t = fabs(uf/f);
  }
  ai = c1 * f - c2 * g;

  return(ai);
}


#undef MAXAIRY
#undef ACC
#undef PI





//AIRY FUNCTION END



// INTEGRAL FUNCTIONS

/************************************************************************
 *                                                                      *
 * Calculates the integral giving the WKB eikonal solution,             *
 * \int^x sqrt(abs(Q(x))) dx                                            *
 * Note chi0 is determined by alpha but is passed to save recomputation *
 * Based on Kosowsky's paper                                            *
 ************************************************************************/
double qintegral_simp(double chi, double chi0, double alpha, int K)
{
  int j, n;
  double integrand, arg, delta, sum, weight;

  /* Simpson's Rule integration for now with n partitions
   */
  n=300;
  delta = (chi - chi0) / (double) n;

  sum = 0.;
  for(j=0; j<= n; j++){
    weight = 1.;
    if(j==0 || j==n) weight = 17./48.;
    if(j==1 || j == n-1) weight = 59./48.;
    if(j==2 || j == n-2) weight = 43./48.;
    if(j==3 || j == n-3) weight = 49./48.;
  
    arg = chi0 + delta * (double)j;
    if(K == -1) integrand = sqrt(fabs(alpha*alpha - 1./(sinh(arg)*sinh(arg))));
    if(K == 0) integrand = sqrt(fabs(alpha*alpha - 1./(arg*arg)));
    if(K == 1) integrand = sqrt(fabs(alpha*alpha - 1./(sin(arg)*sin(arg))));

    sum += weight * integrand;
  }

  sum *= delta;			      

  return(fabs(sum));
}


#define PIOVER2 1.570796327

/************************************************************************
 *                                                                      *
 * Evaluates the exact solution to  the integral giving the WKB         *
 * eikonal solution,   \int^x sqrt(abs(Q(x))) dx                        *
 *                                                                      *
 * In the open case, this integral costs 1 or 2 square roots, an atan   *
 * and a log; its evaluation will be roughly as expensive as the rest   *
 * of the Phi routine. An analytic fit cannot be substantially faster   *
 * because the dependence on alpha of the y-intercept of the linear     *
 * region of the integrand contains a log and an atan, so at best a fit *
 * can only save the computation of the square roots.                   *
 *                                                                      *
 * Note that for the closed case, the variable arg must be between 0    *
 * and alpha; for arg > alpha, symmetry properties reduce the needed    *
 * integral to this case.                                               *
 *                                                                      *
 ************************************************************************/
double qintegral(double sin_K, double alpha, int K)
{
  double exact, arg, root1, root2;

  arg = alpha * sin_K;

  if(K == 0){
    if(arg > 1.){  
      exact = sqrt(arg*arg - 1.) - acos(1./arg);
      return(exact);
    } else {
      root1 = sqrt(1.-arg*arg);
      exact = log((1.+root1)/arg) - root1;
      return(exact);
    }
  }
  else if(K== -1){
    if(arg > 1.){
      root1 = sqrt(arg*arg - 1.);
      root2 = sqrt(arg*arg + alpha*alpha);
      exact = (alpha/2.) * log((2.*arg*arg + alpha*alpha -1. + 2.*root1*root2)
                                 /(1.+alpha*alpha))
                  + atan(root2 / (alpha * root1)) - PIOVER2;
      return(exact);
    } else {  
      root1 = sqrt((1. - arg*arg)*(arg*arg + alpha*alpha));
      exact = (alpha/2.) * atan(-2.*root1/(2.*arg*arg + alpha*alpha -1.))  
              + 0.5*log((4.*alpha*root1 + 4.*alpha*alpha 
                                        + 2.*arg*arg*(1.-alpha*alpha))
                                 /(2.*arg*arg*(1.+alpha*alpha)));
      if(2.*arg*arg + alpha*alpha -1. < 0.) exact -= alpha * PIOVER2;
      return(exact);
    }
  }
  else {
    if(arg > 1.){
      root1 = sqrt(arg*arg - 1.);
      root2 = sqrt(alpha*alpha - arg*arg);
      exact = (alpha/2.) * atan(-2.*root1*root2
                            /(2.*arg*arg - alpha*alpha - 1.))
              - atan(-root2/(alpha*root1))  - PIOVER2;
      if(2.*arg*arg - alpha*alpha -1. > 0.) exact += alpha * PIOVER2;
      return(exact);
    } else {
      root1 = sqrt((1. - arg*arg)*(alpha*alpha - arg*arg));
      exact = atanh(alpha*(1.-arg*arg)/root1)
                - (alpha/2.) * log((alpha*alpha - 2.*arg*arg + 1. + 2.*root1)
                                               /(alpha*alpha - 1.));
      return(exact);
    }
  }			       
}

#undef PIOVER2


// INTEGRAL FUNCTIONS






// HYPERBOLIC BESSEL CLASS BEGIN


#define ACC 40.0
#define BIG 1.e10


/************************************************************************
 *                                                                      *
 * Calculates Phi(l,beta,chi) using recursion on l.                     *
 * See Abbot and Schaefer, ApJ 308, 546 (1986) for needed               *
 * recursion relations and closed-form expressions for l=0,1.           *
 * (Note: Their variable y is the same as chi here.)                    *
 *                                                                      *
 * When the flag direction is negative, downwards recursion on l        *
 * must be used because the upwards direction is unstable to roundoff   *
 * errors. The downwards recursion begins with arbitrary values and     *
 * continues downwards to l=1, where the desired l value is normalized  *
 * using the closed form solution for l=1. (See, e.g., Numerical        *
 * Recipes of Bessel functions for more detail)                         *
 *                                                                      *
 ************************************************************************/
double phi_recursive(int l, int K, double beta, double chi)
{
  int j, direction, lstart;
  double ell, jay, kay, arg, answer;
  double cot_K, sin_K, root_K;
  double phi0, phi1, phi_plus, phi_zero, phi_minus, b_zero, b_minus; 

  ell=(double) l;

  /* Test input values
   */
  if (l < 0) { fprintf(stderr,"Bessel function index l < 0\n"); exit(1); }
  if (beta <= 0.) { fprintf(stderr,"Wavenumber beta < 0\n"); exit(1); }
  if (K < 0) K= -1;
  if (K > 0){
     K=1;
     beta = (double)((int) beta);
     if(beta < 3.){fprintf(stderr,"Wavenumber beta < 3 for K=1\n"); 
                                                                exit(1);}
     if(beta <= ell){fprintf(stderr,"Wavenumber beta <= l for K=1\n"); 
                                                                exit(1);}
   }

  kay = (double) K;
  arg = beta * chi;
  if(K == -1){ 
    cot_K = 1./tanh(chi);
    sin_K = sinh(chi);
    root_K = sqrt(beta*beta + ell*ell);
  }
  if(K == 0){
    cot_K = 1./chi;
    sin_K = chi;
    root_K = beta;
  }
  if(K == 1){
    cot_K = 1./tan(chi);
    sin_K = sin(chi);
    root_K = sqrt(beta*beta - ell*ell);
  }


  /* Closed form solution for l=0
   */
  if(fabs(chi) < 1.e-4) phi0 = 1. - arg*arg/6.;
  else{
    if(arg < 1.e-4) phi0 = chi;
      else  phi0 = sin(arg) / (beta * sin_K);
  }
  if(l==0) return phi0;

  /* Closed form solution for l=1
   */
  if(fabs(chi) < 1.e-4 && (K == -1 || K == 1)){
    if(arg < 1.e-4){
      phi1 = beta * beta * chi / (3. * sqrt(1.+ kay *beta*beta));
    } else {
      phi1 = (sin(arg)/arg - cos(arg))/arg;
    }
  } else if (fabs(arg) < 1.e-4 && K == 0){
    phi1 = arg / 3.;
  } else {
    if (K == -1 || K == 1){
      phi1 = sin(arg) * cot_K / (beta * sin_K) 
                 - cos(arg) / sin_K;
      phi1 /= sqrt(beta * beta - kay);
    }
    if (K == 0) phi1 = (sin(arg)/arg - cos(arg))/arg;
  }
  if(l==1) return phi1;

  /* Find recursion direction
   * direction = +1 for upward recursion, -1 for downward
   */
  if(fabs(cot_K) < root_K / ell) direction = 1;
    else direction = -1;    

if(K==1)direction = 1;

  /* Do upwards recursion on l
   */
  if(direction == 1){
    b_minus = sqrt(beta*beta - kay);
    phi_minus = phi0;
    phi_zero = phi1;

    for(j=2; j<=l; j++){
      jay = (double) j;

      if(K == 0){
        phi_plus = ((2.*jay - 1.) * cot_K*phi_zero - beta*phi_minus)/ beta;
      } else {
        b_zero = sqrt(beta*beta - kay*jay*jay);  
        phi_plus = ((2.*jay - 1.) * cot_K * phi_zero - b_minus * phi_minus) 
                                          / b_zero;
        b_minus = b_zero;
      }
      phi_minus = phi_zero;
      phi_zero = phi_plus;
    }

    return(phi_plus);
  }

  /* Do downwards recursion on l
   */
  else {
    lstart = l + 2 * (int) sqrt(ell*ACC);
    b_zero = sqrt(beta*beta - kay * (double)(lstart*lstart));
    phi_plus = 0.;
    phi_zero = 1.;
    answer = 0.;

    for(j= lstart - 2; j > 0; j--){
      jay = (double) j;

      if(K == 0){
        phi_minus = ((2.*jay + 3.) * cot_K * phi_zero - beta * phi_plus) 
                                           / beta;
      } else {
        b_minus = sqrt(beta*beta - kay*(jay+1.)*(jay+1.));  
        phi_minus = ((2.*jay + 3.) * cot_K * phi_zero - b_zero * phi_plus) 
                                          / b_minus;
        b_zero = b_minus;
      }      
      phi_plus = phi_zero;
      phi_zero = phi_minus;

      if(j == l) answer = phi_minus;
      if(fabs(phi_zero) > BIG){
        phi_plus /= BIG;
        phi_zero /= BIG;
        answer /= BIG;
      }
    }

    /* Normalize answer to previously computed phi1
     */
    answer *= phi1 / phi_minus;
    return(answer);    

  }
}

#undef ACC
#undef BIG

#define PI 3.1415926536
#define ROOTPI 1.772453851
#define ROOT2PI 2.506628275

/************************************************************************
 *                                                                      *
 * Calculates Phi(l,beta,chi) using the Langer uniform approximation    *
 * to the first-order WKB approximation.                                *
 * See C.M. Bender and S.A. Orszag,  Mathematical Methods for           *
 * Scientists and Engineers (McGraw-Hill, 1978; LC QA371.B43),          *
 * chapter 10.                                                          *
 *                                                                      *
 * Differential equation for needed function can be cast into the       *
 * Schrodinger form      \epsilon^2 y'' = Q(x) y                        *
 * where \epsilon^2 = 1/l(l+1) and Q(x) depends on the parameter        *
 * alpha \equiv beta * epsilon.                                         *
 *                                                                      *
 * In the K= +1 case, the function is                                   *
 * determined by its value on the interval [0, pi/2] and the symmetry   *
 * conditions Phi(chi + pi) = (-1)^{beta - l - 1} Phi(chi),             *
 *            Phi(pi - chi) = (-1)^{beta - l - 1} Phi(chi).             *
 * This interval contains one turning point, so the Langer formula      *
 * can be used.                                                         *
 * Note that the second condition at chi = pi/2 gives an eigenvalue     *
 * condition on beta, which  must corrected. For the lowest             *
 * eigenvalue(s), the region between the turning points is not large    *
 * enough for the asymptotic solution to be valid, so the functions     *
 * have a small discontinuity or discontinuous derivative at pi/2;      *
 * this behavior is corrected by employing a 4-term asymptotic          *
 * series around the regular point chi=pi/2.                            *
 * The exact eigenvalue condition requires that beta must be an         *
 * integer >= 3 with beta > l. Note this implies alpha > 1.             *
 *                                                                      *
 ************************************************************************/
double phi_WKB(int l, int K, double beta, double chi)
{
  int betaint;
  double ell, sin_K, symm;
  double transport, eikonal, wkb, arg, arg4;
  double epsilon, alpha, chi0, x, a, b; 
  double asymp, v, v2, v3, A;
  
  ell=(double) l;
  symm = 1.;

  /* Test input values
   */
  if (l < 0) { fprintf(stderr,"Bessel function index l < 0\n"); exit(1); }
  if (beta <= 0.) { fprintf(stderr,"Wavenumber beta < 0\n"); exit(1); }
  if (K < 0) K= -1;
  if (K > 0){
     K=1;
     betaint = (int) beta;
     if(betaint < 3){fprintf(stderr,"Wavenumber beta < 3 for K=1\n"); exit(1);}
     if(betaint <= l){fprintf(stderr,"Wavenumber beta <= l for K=1\n"); 
                                                                    exit(1);}
  }

  /* For closed case, find equivalent chi in [0,pi/2]
   */
  if(K == 1){
    chi = chi - 2.*PI * (int)(chi/(2.*PI));
    if(chi > PI){
      chi = fabs(chi - PI);
      if((betaint-l-1) % 2 == 1) symm *= -1.;
    }
    if(chi > PI/2.){
      chi = PI - chi;
      if((betaint-l-1) % 2 == 1) symm *= -1.;
    }
  }

  /* For closed case, find corrected eigenvalue beta
   */
  if(K== 1) beta = (double) betaint - 1./(8.*ell) + 1./(16.*ell*ell);
   
  /* Definitions
   */
  if(K == -1) sin_K = sinh(chi);
  if(K == 0) sin_K = chi;
  if(K == 1) sin_K = sin(chi);

  epsilon = 1./sqrt(ell*(ell+1.));
  alpha = epsilon * beta;

  /* Calculate the turning point where Q(x)=0.
   * Function in question has only a single simple turning point. 
   */
  if(K == -1) chi0 = asinh(1./alpha);
  if(K == 0) chi0 = 1./alpha;
  if(K == 1) chi0 = asin(1./alpha);

  /* Very close to 0, return 0 (exponentially damped)
   */
  if(fabs(chi) < 1.e-8) return 0.;

  /* Very close to chi0, use usual wkb form to avoid dividing by zero
   */
  if(fabs(chi-chi0) < 1.e-6){

    /* Calculate coefficients of linear and quadratic terms in Q(x) expansion 
     * in the neighborhood of the turning point
     */
    if(K == -1){ 
      a = 2.*alpha*alpha * sqrt(alpha*alpha + 1.);
      b = 3.*alpha*alpha*alpha*alpha + 3.*alpha*alpha;
    }
    if(K == 0){
      a = 2.*alpha*alpha*alpha;
      b = 3.*alpha*alpha*alpha;
    }
    if(K == 1){
      a = 2.*alpha*alpha * sqrt(alpha*alpha - 1.);
      b = 3.*alpha*alpha*alpha*alpha - 3.*alpha*alpha;
    }

    /* Dependent variable x for which Q(x)=0 at x=0
     */
    x = chi0 - chi;

    /* Argument of Airy function
     */
    arg = (x + b*x*x/(5.*a))/pow(epsilon*epsilon/a, 0.3333333333);    

    /* Evaluate Airy function 
     */
    wkb = airy_ai(arg);

    /* Rest of functional dependence
     */
    wkb *= (1.-b*x/(5.*a)) / sin_K;

    /* Normalization factor:
     */
    wkb *= symm * ROOTPI * pow(a*epsilon, -0.166666667) * sqrt(epsilon / beta);

    return(wkb);
  }

      
  /* Langer approximation.
   */
  
  /* Transport factor:
   */
  transport = pow(fabs(1./(sin_K*sin_K) - alpha*alpha), -0.25) / sin_K; 
    
  /* Eikonal factor:
   */
  eikonal = qintegral(sin_K,alpha,K);
  arg = pow(3.*eikonal/(2.*epsilon), 0.1666666667);
  arg4 = arg*arg*arg*arg;
  if(chi > chi0) arg4 *= -1.;

  /* Multiply factors together
   */
  wkb = transport * airy_ai(arg4) * arg;

  /* Normalization factor:
   */
  wkb *= symm * ROOTPI * sqrt(epsilon/beta);

  return(wkb);
}

#undef PI
#undef ROOTPI
#undef ROOT2PI




// HYPERBOLIC BESSEL CLASS END

int main()
{
  int i, l, K;
  double beta, chi, phi1, phi2, jl, jld;
  double epsilon, sqrteps, alpha, chi0, ell, a, sin_K;
 

  FILE *fout;

  fout = fopen("Hyperjl.out","w");

  K = -1;
  
  printf("Enter l and beta");
  scanf("%d %lf",&l,&beta);

  for(i=1; i<=4000; i++){

    ell = (double) l;
    epsilon = 1./sqrt(ell*(ell+1.));
    alpha = epsilon * beta;  

    chi0 = asinh(1./alpha);
    chi = fabs(chi0 - 3. * epsilon / alpha 
                        + (double) i * epsilon/ (20.*alpha));
    if(chi < 0.) continue;
    
    phi1 = phi_recursive(l,K,beta,chi);
    phi2 = phi_WKB(l,K,beta,chi);
    
    printf("%e %f %f\n",chi,phi1,phi2);
    fprintf(fout,"%f %e %e\n",chi,phi1,phi2);
  }


  fclose(fout);
}