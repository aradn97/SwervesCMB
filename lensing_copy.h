#ifndef LENSING_H
#define LENSING_H

#include <stdio.h>
#include <math.h>
//#include "variables.h"
// #include "initval.h"
// #include "others.h"
#include "scalar.h"
#include "numericx.h"

double Sigmal(double);
double C2(double);

class lensing : public initval
{

 public:
#define kmax 5000
#define nknl 50
#define ntheta 50
#define ntaulens 50 // Sets the sampling time for lensing

	double akk[kmax + 1], tf[kmax + 1], tfpr[kmax + 1];
	double atrg[ntaulens + 1][5], tht, rr;
	constexpr static double fourpi = 4 * 3.14159265, pi = 3.14159265;
	double delta2nl[nknl + 1][ntaulens + 1], alnknl[nknl + 1][ntaulens + 1], d2nlp[nknl + 1][ntaulens + 1];
	double alnknlpass[nknl + 1], delta2nlpass[nknl + 1], d2nlppass[nknl + 1] = {0};
	double yfsp1;
	double epsilon[ntheta + 1], theta[ntheta + 1], c2oth[ntheta + 1];

	double fnu, fcb; // (Eisenstein, Hu 1997), Page no 4
	double apcb, yfsok2;
	double xnorm8;
	double omegam;

	int k0;
	int in;

	numericx Numericx;

	/************************************************************************************************************************************/
	/*     rombint returns the integral from a to b of func using Romberg integration. The method converges provided that f(x) is 		*/
	/*     continuous in (a,b). tol indicates the desired relative accuracy in the integral.															*/
	/*		 The function is taken from : http://wise-obs.tau.ac.il/~barkana/InfallSub.c																	*/
	/************************************************************************************************************************************/

	double rombint1(double (*func)(double), double a, double b, double tol)
	{
		int MAXJ = 5;
		int MAXITER = 30, nint, i, j, k, jmax;
		double g[MAXJ + 1], h, gmax, g0, fourj, g1, error;

		h = 0.5 * (b - a);
		gmax = h * (func(a) + func(b));
		g[1] = gmax;
		nint = 1;
		error = 1.e20;
		for (i = 1; !((i > MAXITER) || ((i > 5) && (fabs(error) < tol))); i++)
		{
			// Calculate next trapezoidal rule approximation to integral.
			g0 = 0.0;
			for (k = 1; k <= nint; k++)
			{
				g0 += func(a + (k + k - 1) * h);
				printf("In rombint, %-15.8e\t%-15.8e\n", exp(a + (k + k - 1) * h), g0);
			}

			g0 = 0.5 * g[1] + h * g0;
			h = 0.5 * h;
			nint *= 2;
			jmax = (i < MAXJ) ? i : MAXJ;
			fourj = 1.;
			for (j = 1; j <= jmax; j++)
			// Use Richardson extrapolation.
			{
				fourj *= 4.0;
				g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
				g[j] = g0;
				g0 = g1;
			}
			if (fabs(g0) > tol) //1.e-30)
				error = 1.0 - gmax / g0;
			else
				error = gmax;
			gmax = g0;
			g[jmax + 1] = g0;
		}

		if ((i > MAXITER) && (fabs(error) > tol))
			printf("rombint failed to converge; integral=%g, error=%g\n", g0, error);
		return g0;
	}

	//  This function evaluates the transfer function at wavenumber ak.
	double psl(double alnk)
	{

		double h, t2;
		double apowerspass;
		double ak;
		scalar Scalar;

		h = Variables.h0 / 100.0;

		if (alnk < akk[1])
			t2 = 1.0;
		else
		{
			t2 = Numericx.splint_v2(akk, tf, tfpr, k0, alnk);
			t2 = exp(2.0 * t2);
		}

		ak = exp(alnk) * h;

		Scalar.powersflat(ak, in, &apowerspass);
		return apowerspass * ak / h * t2;
	}

	// Nonlinear mapping from Peacock and Dodds
	// quick fix for WDM, really not appropriate
	// arXiv:astro-ph/9311057v1
	// eqn (23)

	double fnl(double x, double g, double aneff)
	{
		double dum;

		if (aneff > -2.990)
			dum = log(1.0 + aneff / 3.0);
		else
			dum = -5.7;

		double a = 0.482 * exp(-0.947 * dum);
		double b = 0.226 * exp(-1.778 * dum);
		double alpha = 3.310 * exp(-0.244 * dum);
		double beta = 0.862 * exp(-0.287 * dum);
		double v = 11.550 * exp(-0.4230 * dum);
		double dum1 = 1.0 + b * beta * x + exp(alpha * beta * log(a * x));
		double dum2 = 1.0 + exp(beta * (alpha * log(a * x) + 3.0 * log(g) - log(v) - 0.5 * log(x)));
		double fnl1 = x * exp(log(dum1 / dum2) / beta);
		return (fnl1);
	}

	/////////////////
	//     computes integrand of sigma as a function of theta do the time integration at a fixed number of points using preevaluated nonlinear power spectrum;
	//     interpolate to get the ps at a given k

	double sigmal(double alnk)
	{
		double hc = 2.998e3;
		double sigma1;
		double g0 = atrg[ntaulens][4];
		double aux1, aux2, a, r, g;
		double x;
		double yfsp1, ps;
		double aj0, daj0;
		double delta2;
		double dum1[ntaulens];
		double w, dtau;
		double aklens;

		aklens = exp(alnk);
		// EH (97) fitting formula for massive neutrinos
		printf("\nBegin sigmal L165\n");
		if (fnu > 0.0) // arXiv:astro-ph/9710252v1
		{					// eq(12) and (13)
			yfsp1 = yfsok2 * aklens * aklens + 1.0;
			aux1 = exp(0.7 * log(g0 / yfsp1));
			aux1 = aux1 + exp(0.7 * log(fcb) / apcb);
			aux2 = apcb * log(aux1) / 0.7 + (1.0 - apcb) * log(g0);
			g0 = exp(aux2);
		}
		printf("\nIn sigmal L174\n");
		printf("\nIn sigmal ntaulens: %d\n",ntaulens);
		for (int i = 1; i <= ntaulens; i++)
		{
			a = atrg[i][2];
			r = atrg[i][3];
			g = atrg[i][4];
		printf("\nsigmal begin for loop L181 : i, ntaulens : %d %d %lf \n",i,ntaulens,alnk);	

			// EH (97) fitting formula for massive neutrinos
			if (fnu > 0.0) // arXiv:astro-ph/9710252v1
			{					// eq(12) and (13)
				yfsp1 = yfsok2 * aklens * aklens + 1.0;
				aux1 = exp(0.7 * log(g / yfsp1));
				aux1 = aux1 + exp(0.7 * log(fcb) / apcb);
				aux2 = apcb * log(aux1) / 0.7 + (1.0 - apcb) * log(g);
				g = exp(aux2);
			}
			x = aklens * r * tht;
			if (Variables.ilin != 1)
			{
				ps = 9.0 / 4.0 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
			}
			else
			{
				if (alnk < alnknl[1][i])
				{
					ps = 9.0 / 4.0 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
				}
				else //     nonlinear theory next 4 lines
				{
					
					for (int j = 1; j <= nknl; j++)
					{
						alnknlpass[j] = alnknl[j][i];
						delta2nlpass[j] = delta2nl[j][i];
						d2nlppass[j] = d2nlp[j][i];
					}
					printf("\n Line 205 : nknl = %d",nknl);
					delta2 = Numericx.splint_v2(alnknlpass, delta2nlpass, d2nlppass, nknl, alnk);
					ps = 9.0 / 4.0 / fourpi * exp(delta2) * exp(-4.0 * alnk) / hc / hc / hc / hc / a / a * omegam * omegam;
				}
			}
			if (x < 0.1)
			{
				aj0 = 1.0;
				daj0 = x * x / 4.0;
			}
			else
			{
				aj0 = Numericx.bessj0(x);
				daj0 = 1.0 - aj0;
			}

			w = 1.0 - r / rr;
			dum1[i] = fourpi * ps * aklens * w * w * daj0;
		}
		printf("\nIn sigmal L229\n");
		dtau = (Variables.tau0 - Variables.taur) / ntaulens; //     integrate in equal time intervals

		sigma1 = Numericx.integ(dum1, ntaulens);
		sigma1 = sigma1 * dtau;
		printf("\nEnd sigmal L235\n");
		return (sigma1);
	}

	//////////////////////////

	//     computes integrand of sigma as a function of theta do the time integration at a fixed number of points
	//     using preevaluated nonlinear power spectrum; interpolate to get the ps at a given k

	double c2(double alnk)
	{
		//	int ntau=50;
		double hc = 2.998e3;
		double h = Variables.h0 / 100.0;
		double g0 = atrg[ntaulens][4];
		double aux1, aux2, a, r, g;
		double x;
		double yfsp1, ps;
		double delta2;
		double dum1[ntaulens];
		double aklens;
		double aj2, w; //,dtau;

		h = Variables.h0 / 100.0;

		//fourpi=4.0*3.1415926535898;
		aklens = exp(alnk);
		g0 = atrg[ntaulens][4];

		// EH (97) fitting formula for massive neutrinos

		if (fnu > 0.0)
		{
			yfsp1 = yfsok2 * aklens * aklens + 1.0;
			aux1 = exp(0.70 * log(g0 / yfsp1));
			aux1 = aux1 + exp(0.70 * log(fcb) / apcb);
			aux2 = apcb * log(aux1) / 0.70 + (1.00 - apcb) * log(g0);
			g0 = exp(aux2);
		}

		for (int i = 1; i <= ntaulens; i++)
		{
			a = atrg[i][2];
			r = atrg[i][3];
			g = atrg[i][4];

			// EH (97) fitting formula for massive neutrinos
			if (fnu > 0.0)
			{
				yfsp1 = yfsok2 * aklens * aklens + 1.0;
				aux1 = exp(0.7 * log(g / yfsp1));
				aux1 = aux1 + exp(0.70 * log(fcb) / apcb);
				aux2 = apcb * log(aux1) / 0.70 + (1.00 - apcb) * log(g);
				g = exp(aux2);
			}

			x = aklens * r * tht;

			if (Variables.ilin != 1)
				ps = 9.0 / 4.0 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
			else //     linear theory next 2 lines
			{
				if (alnk < alnknl[1][i])
					ps = 9.0 / 4.0 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
				else //     nonlinear theory next 4 lines
				{
					for (int j = 1; j <= nknl; j++)
					{
						alnknlpass[j] = alnknl[j][i];
						delta2nlpass[j] = delta2nl[j][i];
						d2nlppass[j] = d2nlp[j][i];
					}
					delta2 = Numericx.splint_v2(alnknlpass, delta2nlpass, d2nlppass, nknl, alnk);
					ps = 9.0 / 4.0 / fourpi * exp(delta2) * exp(-4.0 * alnk) / hc / hc / hc / hc / a / a * omegam * omegam;
				}
			}

			if (x < 0.1)
				aj2 = x * x / 8.0;
			else
				aj2 = Numericx.bessj(2, x);

			w = 1.0 - r / rr;
			dum1[i] = fourpi * ps * aklens * w * w * aj2;
		}
		double dtau = (Variables.tau0 - Variables.taur) / ntaulens; //     integrate in equal time intervals
		double c2n = Numericx.integ(dum1, ntaulens);
		c2n = c2n * dtau;
		return (c2n);
	}

	//     computes integrand of sigma as a function of theta do the time integration at a fixed number of points
	//     using preevaluated nonlinear power spectrum; interpolate to get the ps at a given k

	double psk(double xkappa)
	{
		double hc = 2.998e3;
		double h = Variables.h0 / 100.0;
		double g0 = atrg[ntaulens][4];
		double aux1, aux2, a, r, g;
		double x;
		double yfsp1, ps;
		double delta2;
		double dum1[ntaulens];
		double alnk, aklens;
		double w;

		h = Variables.h0 / 100.0;

		g0 = atrg[ntaulens][4];

		// EH (97) fitting formula for massive neutrinos

		if (fnu > 0.0)
		{
			yfsp1 = yfsok2 * aklens * aklens + 1.00;
			aux1 = exp(0.70 * log(g0 / yfsp1));
			aux1 = aux1 + exp(0.70 * log(fcb) / apcb);
			aux2 = apcb * log(aux1) / 0.70 + (1.0 - apcb) * log(g0);
			g0 = exp(aux2);
		}

		for (int i = 1; i <= ntaulens - 1; i++)
		{
			a = atrg[i][2];
			r = atrg[i][3];
			g = atrg[i][4];
			aklens = xkappa / r;
			alnk = log(aklens);

			// EH (97) fitting formula for massive neutrinos
			if (fnu > 0.0)
			{
				yfsp1 = yfsok2 * aklens * aklens + 1.0;
				aux1 = exp(0.70 * log(g / yfsp1));
				aux1 = aux1 + exp(0.70 * log(fcb) / apcb);
				aux2 = apcb * log(aux1) / 0.7 + (1.0 - apcb) * log(g);
				g = exp(aux2);
			}

			x = aklens * r * tht;

			if (Variables.ilin != 1)
			{
				ps = 9.00 / 4.00 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
			}
			else
			{ // linear theory next 2 lines
				if ((alnk < alnknl[1][i]) || (alnk > alnknl[nknl][i]))
				{
					ps = 9.0 / 4.00 * xnorm8 * psl(alnk) / exp(alnk) * g * g / g0 / g0 / hc / hc / hc / hc / a / a * omegam * omegam;
				}
				else // nonlinear theory next 4 lines
				{
					for (int j = 1; j <= nknl; j++)
					{
						alnknlpass[j] = alnknl[j][i];
						delta2nlpass[j] = delta2nl[j][i];
						d2nlppass[j] = d2nlp[j][i];
					}
					delta2 = Numericx.splint_v2(alnknlpass, delta2nlpass, d2nlppass, nknl, alnk);
					ps = 9.00 / 4.00 / fourpi * exp(delta2) * exp(-4.00 * alnk) / hc / hc / hc / hc / a / a * omegam * omegam;
				}
			}

			w = 1.0 - r / rr;

			dum1[i] = fourpi * ps * aklens * w * w;
		}

		dum1[ntaulens] = 0;

		double dtau = (Variables.tau0 - Variables.taur) / ntaulens; //   integrate in equal time intervals//
		double pskn = Numericx.integ(dum1, ntaulens);
		pskn = pskn * dtau;
		return (pskn);
	}

#define nx 7 * l0max + 2

	//  calculating the multipoles C_l after lensing using sigma(theta)

	void cllens(int in)
	{
		double clt[l0max], cle[l0max], clb[l0max], clc[l0max];
		double cltl[l0max] = {0.0}, clel[l0max] = {0.0}, clbl[l0max] = {0.0}, clcl[l0max] = {0.0};
		//	double clkk[l0max],cltk[l0max];

		int nth = 50;

		double aj0[nx], daj0[nx];
		double aj2[nx], daj2[nx];
		double aj4[nx], daj4[nx];
		double aj6[nx], daj6[nx];

		double c2othprime[nth + 1], epsprime[nth];

		double epslo, epshi, dx;
		int ith, nxfile, nt;
		int lmaximfile;
		double xmaxlens;

		int l0, m;
		double dtheta;
		double sig1, c2gl, cth1, cth2, cth3, cth4;
		double theta1, eps1;
		double x1, b, b0, b2, b4, b6;
		double aux1, aux2;
		double a31, b31;
		int i1lo, i1hi;
		double x2;

		int lmaxin = Variables.lmo;

		numericx Numericx;

		FILE *fwrite;
		FILE *fplens;
		fwrite = fopen("cllens.txt", "w");

		// Getting Cls
		for (int l = 2; l <= lmaxin; l++)
		{
			clt[l] = GlobalArray.clts[l][in] / l / (l + 1.0);
			cle[l] = GlobalArray.cles[l][in] / l / (l + 1.0);
			clc[l] = GlobalArray.clcs[l][in] / l / (l + 1.0);
		}
		// Get interpolation tables for epsilon and c2.

		ith = ntheta;
		theta[ith] = pi;
		epsilon[ith] = 0.0;
		c2oth[ith] = 0.0;
		epslo = 1.e40;
		epshi = 1.e40;

		Numericx.spline(theta, epsilon, ith, epslo, epshi, epsprime);
		Numericx.spline(theta, c2oth, ith, epslo, epshi, c2othprime);

		//     Read array with bessel function

		fplens = fopen("jlens.dat", "r");

		fscanf(fplens, "%d %d %d %lf", &nxfile, &lmaximfile, &nt, &dx);

		for (int i = 1; i <= nxfile; i++)
		{
			fscanf(fplens, "%lf %lf %lf %lf %lf %lf %lf %lf", &aj0[i], &daj0[i], &aj2[i], &daj2[i], &aj4[i], &daj4[i], &aj6[i], &daj6[i]);
		}

		if (nxfile != nx)
		{
			printf("\njlens.dat and lensing.h have inconsistent");
			printf("\nl0max, you will have to rerun jlens.dat");
		}
		xmaxlens = pi * lmaxin / 10.0;

		if (lmaxin > lmaximfile)
		{
			printf("\njlens.dat was was run with a lower lmax");
			printf("\nthan needed. You will have to rerun it");
		}

		//     No need to integrate to pi because we are calculating the deviations due to lensing which are very small on large
		//     angular scales

		l0 = 2 * lmaxin;
		m = l0 / 10;
		dtheta = pi / l0;

		for (int i = 1; i <= m; i++)
		{
			theta1 = i * dtheta;

			eps1 = Numericx.splint_v2(theta, epsilon, epsprime, ith, theta1);
			sig1 = pow((eps1 * theta1), 2);

			c2gl = Numericx.splint_v2(theta, c2oth, c2othprime, ith, theta1);
			c2gl = c2gl * theta1 * theta1;

			//     get the correlation functions at theta
			cth1 = 0.0;
			cth2 = 0.0;
			cth3 = 0.0;
			cth4 = 0.0;

			for (int l1 = 2; l1 <= lmaxin; l1++)
			{
				x1 = l1 * theta1; //  Interpolate Bessel function

				i1lo = (int)(x1 / dx) + 1;
				i1hi = i1lo + 1;
				a = i1lo - 1.0 - x1 / dx;
				b = 1.0 - a;
				a31 = (a * a * a - a) * dx * dx / 6.0;
				b31 = (b * b * b - b) * dx * dx / 6.0;

				b0 = a * aj0[i1lo] + b * aj0[i1hi] + a31 * daj0[i1lo] + b31 * daj0[i1hi];
				b2 = a * aj2[i1lo] + b * aj2[i1hi] + a31 * daj2[i1lo] + b31 * daj2[i1hi];
				b4 = a * aj4[i1lo] + b * aj4[i1hi] + a31 * daj4[i1lo] + b31 * daj4[i1hi];
				b6 = a * aj6[i1lo] + b * aj6[i1hi] + a31 * daj6[i1lo] + b31 * daj6[i1hi];

				aux1 = exp(-0.5 * sig1 * (l1 * l1));
				aux2 = 0.5 * (l1 * l1) * c2gl * aux1;
				aux1 = x1 * (aux1 - 1.0);
				aux2 = x1 * aux2;

				//     cth4 is - the correlation function, when we go back to fourier space that sign is compensated.

				cth1 = cth1 + (aux1 * b0 + aux2 * b2) * clt[l1];
				cth2 = cth2 + (aux1 * b0 + aux2 * b2) * cle[l1];
				cth3 = cth3 + (aux1 * b4 + aux2 * 0.5 * (b2 + b6)) * cle[l1];
				cth4 = cth4 + (aux1 * b2 + aux2 * 0.5 * (b0 + b4)) * clc[l1];
			}

			for (int l = 2; l <= lmaxin; l++)
			{
				x2 = l * theta1;

				i1lo = (int)(x2 / dx) + 1;
				i1hi = i1lo + 1;
				a = i1lo - 1.0 - x2 / dx;
				b = 1.0 - a;
				a31 = (a * a * a - a) * dx * dx / 6.0;
				b31 = (b * b * b - b) * dx * dx / 6.0;

				b0 = a * aj0[i1lo] + b * aj0[i1hi] + a31 * daj0[i1lo] + b31 * daj0[i1hi];
				b2 = a * aj2[i1lo] + b * aj2[i1hi] + a31 * daj2[i1lo] + b31 * daj2[i1hi];
				b4 = a * aj4[i1lo] + b * aj4[i1hi] + a31 * daj4[i1lo] + b31 * daj4[i1hi];

				cltl[l] = cltl[l] + b0 * cth1;
				clel[l] = clel[l] + 0.5 * (b0 * cth2 + b4 * cth3);
				clbl[l] = clbl[l] + 0.5 * (b0 * cth2 - b4 * cth3);
				clcl[l] = clcl[l] + b2 * cth4;
			}
		}

		//     Last 300 ls are not accurate because of convolution.

		for (int l = 2; l <= lmaxin - 300; l++)
		{
			aux1 = ((l + 1) * l);
			aux2 = aux1 * dtheta;

			clt[l] = cltl[l] * aux2 + clt[l] * aux1;
			cle[l] = clel[l] * aux2 + cle[l] * aux1;
			clb[l] = clbl[l] * aux2;
			clc[l] = clcl[l] * aux2 + clc[l] * aux1;

			fprintf(fwrite, "%d %e %e %e %e %e %e\n", l, clt[l], cle[l], clb[l], clc[l], GlobalArray.clkk[l][in], GlobalArray.cltk[l][in]);
		}
		fclose(fplens);
	}
#undef nx
	//#undef ntaulens
	//     fluctuations in angular separation due to the potential fluctuations
	//     from Limber's eqn. using nonlinear evolution. normalize to sigma8
	//     units are Mpc

#define nx 1000
	//#define ntaulens 50					 // Sets the sampling time for lensing

	void epsilongen()
	{
		double h;		// Houbble constant in 100km/sec/MpC
		double zeqp1;  // Redshift at the matter and radiation equality
		double apcb;	// (Eisenstein, Hu 1997), Page no 4, eq (11). 'pcb' is defined as 'apcb'
		double yfsok2; // Free streaming epoch as a function of scale. (Eisenstein, Hu 1997), Page no 6, eq (14).
		double akkhi, akklo;
		double aktoq;
		double tcdm, tfb, tfg, tfn, tfnm;
		double a0, a1;
		double alnkmin, alnkmax;
		double ar; //taur;

		//xapl0,
		double xaplo, xtau[nx + 1], xa[nx + 1], xaphi;
		double xpr[nx + 1] = {0};
		double tau, r, x, f, g;
		double omnow, oqnow, term;

		double alnklo, alnkhi, delta2crit, alnk, delta2, alnkcrit, dlnk, g0, g0c;
		double gc, aux1, aux2, aux3, alnklow; //,alnk1,alnlkhigh;
		double d2nlplo, d2nlphi, alnkhigh, aneff, alnkl;
		double t2high, t2low;
		double theta0, theta1, dlntheta;
		double sigth, xkappa, c2th; //,sigma2kappa;
		double akkt, ak;

		double delta = 1.0;

		int i;

		numericx Numericx;

		FILE *fptrans;
		fptrans = fopen("transfcmb.d", "r");
		rewind(fptrans);

		omegam = Variables.omegac + Variables.omegab + Variables.omegan;
		h = Variables.h0 / 100.0;

		zeqp1 = (2.5e4) * omegam * h * h * pow((2.7 / Variables.tcmb), 4); // Redshift at the matter and radiation equality
																								 // arXiv:astro-ph/9710252v1 (Eisenstein, Hu 1997)
		fnu = Variables.omegan / omegam;
		fcb = (Variables.omegac + Variables.omegab) / omegam;

		if (fnu > 0.0)
		{
			apcb = 0.25 * (5.0 - sqrt(1.0 + 24.0 * fcb)); // (Eisenstein, Hu 1997), Page no 4, eq (11). 'pcb' is defined as 'apcb'
			aktoq = (2.70 / Variables.tcmb) * (2.70 / Variables.tcmb) / (omegam * h * h) * h;
			yfsok2 = 17.2 * fnu * (1.0 + 0.488 * pow(fnu, -7.0 / 6.0)) * pow((Variables.annunr * aktoq / fnu), 2); // Free streaming epoch as a function of scale. (Eisenstein, Hu 1997), Page no 6, eq (14).
		}
		else
		{
			apcb = 1.0;
			yfsok2 = 0.0;
		}

		//     We will use transfer function output from CMBFAST and fitting formulae for the growth factor. If more accuracy is needed
		//     one can obtain the transfer function as a function of time numerically from CMBFAST.
		i = 1;

		printf("\nDone done\n"); //exit(1);
		if (Variables.omegan > 0.0)
		{
			for (i = 1; akk[i] < 1.0 && fscanf(fptrans, "%lf %lf %lf %lf %lf %lf", &akkt, &tcdm, &tfb, &tfg, &tfn, &tfnm) != EOF; i++) //akk[i]<=5.0
			{
				akk[i] = log(akkt);
				tf[i] = (Variables.omegac * tcdm + Variables.omegab * tfb + Variables.omegan * tfnm) / omegam;
				tf[i] = log(tf[i]);
			}
		}
		else
		{
			for (i = 1; akk[i] <= 5.0 && fscanf(fptrans, "%lf %lf %lf %lf %lf", &akkt, &tcdm, &tfb, &tfg, &tfn) != EOF; i++) //akk[i]<=5.0
			{
				akk[i] = log(akkt);
				tf[i] = (Variables.omegac * tcdm + Variables.omegab * tfb) / omegam;
				tf[i] = log(tf[i]);
			}
		}
		printf("\nEnd");
		fclose(fptrans);
		k0 = i - 1;

		akkhi = 1.0e40;
		akklo = 1.0e40;

		Numericx.spline(akk, tf, k0, akkhi, akklo, tfpr); // Check what this function is doing

		// Conformal time today
		a0 = 1.0e-8;
		a1 = 1.0;

		Variables.tau0 = h * Variables.tau0;
		ar = 1.0 / 1100.0;
		Variables.taur = h * Numericx.rombint(Dtauda, a0, ar, 1.0e-5);

		alnkmin = akk[1];
		alnkmax = akk[k0];

		// Set up (a,tau,r,g,n_eff) table for ntime points equally spaced in time; use interpolation between a and tau
		xa[0] = 1.0e-8;
		xtau[0] = 0.0;
		for (int i = 1; i <= nx; i++) //nx
		{
			xa[i] = (1.0 * i) / (nx); // nx is defined. So logical substitution.
			xtau[i] = xtau[i - 1] + h * Numericx.rombint(Dtauda, xa[i - 1], xa[i], 1.0e-5);
		}

		xaplo = 1.0 / Dtauda(xa[1]);
		xaphi = 1.0 / Dtauda(1.0);

		Numericx.spline(xtau, xa, nx, xaplo, xaphi, xpr);

		for (int i = 1; i <= ntaulens; i++)
		{
			tau = i * (Variables.tau0 - Variables.taur) / ntaulens + Variables.taur;
			a = Numericx.splint_v2(xtau, xa, xpr, nx, tau);
			r = Variables.tau0 - tau;

			if (r < 0.0)
				r = 0.0;
			x = 1.0 + omegam * (1.0 / a - 1) + Variables.omegav * (a * a - 1.0);
			f = exp(4.0 * log(omegam / a / x) / 7.0);

			//     Fitting formula for linear growth of perturbations. Extra (1+zeq)=zeqp1 so that same normalization as EH 97 to use
			//     their fitting formula for massive neutrino growth rate cancels in other cases.

			g = zeqp1 * 2.5 * omegam / (x * f + 1.5 * omegam / a + 1.0 - omegam - Variables.omegav);

			//     Ma, Caldwell, Bode, Wang fitting formula

			if (Variables.wdyn < 0.0)
			{
				omnow = omegam / (omegam + Variables.omegav * pow(a, (-3.0 * Variables.wdyn)));
				oqnow = 1.0 - omnow;

				term = -(0.255 * oqnow + 0.366 * log(omnow)) - (0.00266 * oqnow - 0.07 * log(omnow)) / Variables.wdyn - (0.305 * oqnow + 0.266 * log(omnow)) * Variables.wdyn;
				g = g * pow((-Variables.wdyn), term); //   arXiv:astro-ph/9906174v1
			}													  //   Ma, Caldwellb, Bodec, Wang (1999)

			atrg[i][1] = tau;
			atrg[i][2] = a;
			atrg[i][3] = r;
			atrg[i][4] = g;
		}

		rr = Variables.tau0 - Variables.taur;

		//    If non linear calculation requested, nonlinear mapping of power spectrum

		printf("\nVariables.ilin = %d", Variables.ilin);
		if (Variables.ilin == 1)
		{
			//    only k>kcrit for which linear Delta^2>0.1  today.
			//    find kcrit between kmin and kmax.
			alnklo = alnkmin;
			alnkhi = alnkmax;
			delta2crit = 0.1;

			while ((alnkhi - alnklo) > 0.01)
			{
				alnk = (alnkhi + alnklo) / 2.0;
				delta2 = fourpi * xnorm8 * exp(alnk * 3.0) * psl(alnk);

				if (delta2 > delta2crit)
					alnkhi = alnk;
				else
					alnklo = alnk;
			}

			alnkcrit = alnk;
			dlnk = (alnkmax - alnkcrit) / (nknl - 1.0);

			g0 = atrg[ntaulens][4];
			g0c = g0;

			printf("\nL770\n");
			for (int i = 1; i <= ntaulens; i++)
			{
				g = atrg[i][4];
				gc = g;
				a = atrg[i][2];

				for (int j = 1; j <= nknl; j++)
				{
					alnkl = alnkcrit + (j - 1.0) * dlnk; // EH (97) fitting formula for massive neutrinos
					if (fnu > 0.0)								 // arXiv:astro-ph/9710252v1
					{												 // Eisenstein, Hu (1997) eq(11) - eq(13)
						ak = exp(alnkl);
						yfsp1 = yfsok2 * ak * ak + 1.00;
						aux1 = exp(0.70 * log(g / yfsp1));
						aux3 = exp(0.70 * log(fcb) / apcb);
						aux1 = aux1 + aux3;
						aux2 = apcb * log(aux1) / 0.70 + (1.00 - apcb) * log(g);
						gc = exp(aux2);

						aux1 = exp(0.70 * log(g0 / yfsp1));
						aux1 = aux1 + aux3;
						aux2 = apcb * log(aux1) / 0.7 + (1.0 - apcb) * log(g0);
						g0c = exp(aux2);
					}
					delta2 = fourpi * xnorm8 * exp(3 * alnkl) * psl(alnkl) * gc * gc / g0c / g0c;
					alnklow = alnkl - delta - log(2.0);
					t2low = log(psl(alnklow));
					alnkhigh = alnkl + delta - log(2.0);
					t2high = log(psl(alnkhigh));
					aneff = (t2high - t2low) / 2.0 / delta;
					delta2nl[j][i] = fnl(delta2, g / a / zeqp1, aneff);
					alnknl[j][i] = log(1.0 + delta2nl[j][i]) / 3.0 + alnkl;

					//     use log(delta2) to facilitate interpolation
					delta2nl[j][i] = log(delta2nl[j][i]);
				}
				d2nlplo = 1.0e40;
				d2nlphi = 1.0e40;

				for (int j = 1; j <= nknl; j++)
				{
					alnknlpass[j] = alnknl[j][i];
					delta2nlpass[j] = delta2nl[j][i];
				}
				Numericx.spline(alnknlpass, delta2nlpass, nknl, d2nlplo, d2nlphi, d2nlppass);
				for (int j = 1; j <= nknl; j++)
					d2nlp[j][i] = d2nlppass[j];
			}
		}

		//    k integration

		printf("\nL823\n");
		theta0 = pi / (180.0 * 600.0);
		theta1 = 0.80;
		dlntheta = (log(theta1) - log(theta0)) / (ntheta - 1);
		printf("\nL827\n");
		printf("\nntheta: %d\n", ntheta);
		for (int itheta = 1; itheta <= ntheta; itheta++) //ntheta
		{
			printf("\nitheta: %d\n", itheta);
			theta[itheta] = exp(log(theta0) + dlntheta * (itheta - 1));
			tht = theta[itheta];
			printf("\ntht: %d\n", tht);
			//sigth=sqrt(fourpi*rombint1(Sigmal,alnkmin,alnkmax,1.0e-5));
			Sigmal(alnkmax);
		/*	printf("\nL837 alnkmin,alnkmax %lf %lf\n",alnkmin, alnkmax);	
			exit(1);		
			sigth = sqrt(fourpi * Numericx.rombint(Sigmal, alnkmin, alnkmax, 1.0e-5));
			printf("\nsigth: %d\n", sigth);
			xkappa = 1 / tht;
			printf("\nxkappa: %d\n", xkappa);
			epsilon[itheta] = sigth / theta[itheta];
			//c2th=(fourpi*rombint1(C2,alnkmin,alnkmax,1.0e-5));
			c2th = (fourpi * Numericx.rombint(C2, alnkmin, alnkmax, 1.0e-5));
			c2oth[itheta] = c2th / theta[itheta] / theta[itheta];
			*/
		}
		printf("\nL86\n");
	}

#undef nx

	void CalcLensing()
	{ // Decide if linear or non-linear
		double xn;
		Variables.ilin = 0;
		if (Variables.lensflag == 2)
			Variables.ilin = 1; // Loop over power spectra
		for (in = 1; in <= Variables.nn; in++)
		{
			xn = Variables.an[in];								// Take the last transfer function as the one for today.
			xnorm8 = Variables.anorm[Variables.ntf][in]; // Calculate epsilon and c2
			printf("%e %e - going epsilon", xn, xnorm8);
			epsilongen(); // Lens the Cls
			cllens(in);
		}
	}
};

double Sigmal(double alnk)
{
	lensing Lensing;
	return Lensing.sigmal(alnk);
}

double C2(double alnk)
{
	lensing Lensing;
	return Lensing.c2(alnk);
}
#endif
