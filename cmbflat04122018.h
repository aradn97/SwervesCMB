/*
 * cmbflat.h
 *
 *  Created on: Mar 5, 2018
 */

#ifndef CMBFLAT_H_
#define CMBFLAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>

#include "initval.h"
#include "variables.h"
#include "numericx.h"
//#include "neutrino.h"
#include "others.h"
#include "reion.h"
#include "recombination.h"
#include "scalar.h"
#include "tensor.h"

class cmbflat : public initval
{
  private:
    /******************************************************************************************************************/
    /*		It will printout the transfer function at the given points.																	*/
    /******************************************************************************************************************/

    double apowers;
    int kmaxjl;

    double **ajl;   //[ketamax0 + 1][lmax + 1]
    double **ajlpr; //[ketamax0 + 1][lmax + 1]	BESSEL'S Function
    double xx[ketamax0 + 1], dxx[ketamax0 + 1];

  public:
    double **matrix(int col, int row)
    {
        double **m;
        m = (double **)malloc(col * sizeof(double));
        if (!m)
        {
            printf("allocation failure 1 in dmatrix()");
            exit(1);
        }

        /* allocate rows and set pointers to them */
        for (int i = 0; i <= col; i++)
        {
            m[i] = (double *)malloc(row * sizeof(double));
            if (!m[i])
            {
                printf("allocation failure 2 in dmatrix()");
                exit(1);
            }
        }
        //printf("Successfully allocated\n");
        /* return pointer to array of pointers to rows */
        return m;
    }

    /*************************************************************************************************************/
    /* 	This function will sort the k arrays in an increasing order														 */
    /*************************************************************************************************************/

    int indexx1(int n, double arr[], int indx[])
    {
        int jstack, l, ir, indxt, i, itemp, k;
        double asan;
        int M = 7, NSTACK = 50;
        int istack[NSTACK + 1];
        int j;
        //	int i;
        for (int js = 1; j <= n; j++)
            indx[js] = js;

        jstack = 0;
        l = 1;
        ir = n;
        //   printf("%d\n",n);
        //   for(int ii=1;ii<=n;ii++)
        //  		printf("%e %d\n",arr[ii],indx[ii]);
        //  exit(1);

        while (1)
        {
            if (ir - l < M)
            {
                //		printf("ir,l = %d %d",ir,l);
                for (int j = l + 1; j <= ir; j++)
                {
                    indxt = indx[j];
                    asan = arr[indxt];
                    i = j - 1;
                    for (i = j - 1; i >= 1; i--)
                    {
                        if (arr[indx[i]] <= asan)
                            break;
                        indx[i + 1] = indx[i];
                        //            printf("%d %d\n",i,j-1);
                    }
                    //         printf("Hi %d",i);
                    if (arr[indx[i]] > asan)
                        i = 0;
                    //         printf("Cought");
                    indx[i + 1] = indxt;
                }

                //		printf("ir,l = %d %d",ir,l);

                if (jstack == 0)
                    return 0;

                ir = istack[jstack];
                l = istack[jstack - 1];
                jstack = jstack - 2;
                // 		printf("%d %d %d\n",ir,l,jstack);
            }

            else
            {
                k = (l + ir) / 2;
                itemp = indx[k];
                indx[k] = indx[l + 1];
                indx[l + 1] = itemp;

                if (arr[indx[l + 1]] > arr[indx[ir]])
                {
                    itemp = indx[l + 1];
                    indx[l + 1] = indx[ir];
                    indx[ir] = itemp;
                }

                if (arr[indx[l]] > arr[indx[ir]])
                {
                    itemp = indx[l];
                    indx[l] = indx[ir];
                    indx[ir] = itemp;
                }

                if (arr[indx[l + 1]] > arr[indx[l]])
                {
                    itemp = indx[l + 1];
                    indx[l + 1] = indx[l];
                    indx[l] = itemp;
                }

                i = l + 1;
                j = ir;
                indxt = indx[l];
                asan = arr[indxt];

            SAN3:
                i = i + 1;

                if (arr[indx[i]] < asan)
                    goto SAN3;

                do
                {
                    j = j - 1;
                } while (arr[indx[j]] > asan);

                if (j >= i)
                {

                    itemp = indx[i];
                    indx[i] = indx[j];
                    indx[j] = itemp;
                    goto SAN3;
                }
                indx[l] = indx[j];
                indx[j] = indxt;
                jstack = jstack + 2;

                if (jstack > NSTACK)
                {
                    printf("NSTACK too small in indexx");
                    exit(1);
                }

                if (ir - i + 1 >= j - l)
                {
                    istack[jstack] = ir;
                    istack[jstack - 1] = i;
                    ir = j - 1;
                }
                else
                {
                    istack[jstack] = j - 1;
                    istack[jstack - 1] = l;
                    l = i;
                }
            }
        }
        return 0;
    }

    /*************************************************************************************************************/
    /*     This subroutine reads the jl files from disk and initializes other variables needed in CMBANS.		 */
    /*************************************************************************************************************/

    void initjl()
    {
        int kmaxfile, lmofile, l0file;
        int ketamax = 3 * l0max + 126, lfile[lmax + 1];
        double xlimmin = 10.0;
        double x, xlim;
        double d0hi = 1.0e40, d0lo = 1.0e40;
        double aj = 0.0;
        FILE *fp;
        numericx Numericx;
        ajl = matrix(ketamax0 + 1, lmax + 1);
        ajlpr = matrix(ketamax0 + 1, lmax + 1);

        fp = fopen("jlgen.dat", "r");

        fscanf(fp, "%d %d", &lmofile, &kmaxfile);
        // The code needs 300 more l's for the lensing calculation.
        if ((Variables.lmo > lmofile) || (Variables.akmax0 > kmaxfile))
        {
            printf("\n\nThe value of lmax and/or kmax is more than what is stored in the file");
            printf("\nIn the file we have lmax=%d and kmax = %d", lmofile - 300, kmaxfile);
            printf("\nSo restart the calculation please. \n");
            exit(1);
        }

        kmaxfile = kmaxfile - 25 + 151;
        kmaxjl = (Variables.akmax0 - 25 + 151);

        double bjl[kmaxfile + 1], bjlpr[kmaxfile + 1];

        if (kmaxfile > ketamax)
        {
            printf("kmax in file %d is too large. Maximum value for ketamax can be %d", kmaxfile, ketamax);
            exit(1);
        }

        fscanf(fp, "%d", &l0file); // Checking if the lvalues.inc file used to build jl file is the same as the one in the code.
        for (int j = 1; j <= Variables.l0; j++)
        {
            fscanf(fp, "%d", &lfile[j]);
            //printf("%d\n", lfile[j]);
        }

        for (int j = 1; j <= Variables.l0; j++)
        {
            if (Variables.l[j] != lfile[j])
            {
                printf("\n\nlvalues.inc file used to build jl file and the one in the code differ.");
                printf("\nYou must use the same one %d %d %d", j, Variables.l[j], lfile[j]);
                exit(1);
            }
        }

        for (int i = 1; i <= kmaxfile; i++) // reading  j_l. remember to create jl.dat with jlgen.f first using correct lmax and ketamax
        {
            if (i < 151)
            {
                if (i <= 51)
                    xx[i] = double(i - 1) / 10.0;
                else
                    xx[i] = double(i - 51) / 5.0 + 5.0;
            }
            else
                xx[i] = (i - 151) + 25.0;
        }

        for (int j = 1; j <= l0file; j++)
        {
            for (int i = 1; i <= kmaxfile; i++)
            {
                x = xx[i];
                xlim = (Variables.l[j]) / 20.0;
                xlim = Numericx.max(xlim, xlimmin);

                xlim = Variables.l[j] - xlim;
                if (x - xlim > 1e-10)
                {
                    aj = 0.0;
                    fscanf(fp, "%lf", &aj);
                    ajl[i][j] = aj;
                }
                else
                    ajl[i][j] = 0.0;
            }
        }
        fclose(fp);

        for (int i = 2; i <= (kmaxfile - 1); i++)
        {
            dxx[i] = (xx[i + 1] - xx[i - 1]) / 2.0;
        }

        dxx[1] = xx[2] / 2.0;
        dxx[kmaxfile] = (xx[kmaxfile] - xx[kmaxfile - 1]) / 2.0;

        bjl[0] = 0;
        bjlpr[0] = 0;

        for (int j = 1; j <= Variables.l0; j++) // Get the interpolation matrix for bessel functions
        {
            for (int k = 1; k <= kmaxfile; k++)
            {
                bjl[k] = ajl[k][j];
            }
            Numericx.spline(xx, bjl, kmaxfile, d0lo, d0hi, bjlpr);

            for (int k = 1; k <= kmaxfile; k++)
                ajlpr[k][j] = bjlpr[k];
        }
    }

    void output_power(int ntf, double amnu)
    {
        int kmax1 = 1000;
        int nmax = 1000;
        int indx[1001] = {0};
        numericx Numericx;

        double k[kmax1 + 1], tfc[kmax1 + 1], tfb[kmax1 + 1], tfg[kmax1 + 1], tfn[kmax1 + 1], tfnm[kmax1 + 1];
        FILE *fp;
        fp = fopen("transfcmb.d", "w");
        for (int itf = 1; itf <= ntf; itf++)
        {
            if (GlobalArray.istore[itf] > nmax)
            {
                printf("\n Need to increase nmax");
                exit(1);
            }
            for (int ncount = 1; ncount <= GlobalArray.istore[itf]; ncount++)
            {
                k[ncount] = GlobalArray.power[6 * (ncount - 1) + 1][itf];
                tfc[ncount] = GlobalArray.power[6 * (ncount - 1) + 2][itf];
                tfb[ncount] = GlobalArray.power[6 * (ncount - 1) + 3][itf];
                tfg[ncount] = GlobalArray.power[6 * (ncount - 1) + 4][itf];
                tfn[ncount] = GlobalArray.power[6 * (ncount - 1) + 5][itf];
                tfnm[ncount] = GlobalArray.power[6 * (ncount - 1) + 6][itf];
            }

            /*************************************************************************************************************/
            /*  	sort k's in increasing order																									 */
            /*************************************************************************************************************/

            Numericx.indexx(GlobalArray.istore[itf], k, indx);
            for (int icount = 1; icount <= GlobalArray.istore[itf]; icount++)
            {
                if (k[indx[icount]] != 0.0)
                {
                    if (amnu != 0.0)
                        fprintf(fp, "%e %e %e %e %e %e\n", k[indx[icount]], tfc[indx[icount]], tfb[indx[icount]], tfg[indx[icount]], tfn[indx[icount]], tfnm[indx[icount]]);
                    else
                        fprintf(fp, "%e %e %e %e %e\n", k[indx[icount]], tfc[indx[icount]], tfb[indx[icount]], tfg[indx[icount]], tfn[indx[icount]]);
                }
            }
        }
        fclose(fp);
    }

    void CMBflat()
    {
        numericx Numericx;
        neutrino Neutrino;
        reion Reion;
        recombination Recombination;

        double omegam; // Matter density
        double zst;    // Redshift upto which the equations will be integrated

        double arel = 0.0; //mc^2/kT for neutrinos

        double akmax, akmin, dlnk; // akmin -> k from which the integration will start
                                   // akmax -> integration will be done till this k value
                                   // dlnk is the logarithmic spacing used for low k. Check P6

        double *ak0; //Store the wavenumbers at the gridpoints

        double taumax, taumin, dlntau0, taurend; // taumax -> where to stop integration
                                                 // taumin -> where to start integration
                                                 // taurend -> Recombination End time
                                                 // dlntau0 -> log time step

        double shor; //sound horizon

        double a0; // Starting scale factor

        double aprevious;
        double *tautf0;                      // Temp array used for calculating time steps
        double tautf[ntau + ntfmax + 1];     // tautf <- stores the conformal time at which the transfer functions are requested + the conformal time where we need the transfer functions for lensing
        int indtau[ntau + ntfmax + 1] = {0}; // It will store the indices of the points where the transfer functions are to be calculated

        double tauend, taulasttf, taustart, taustart0;
        double tau;

        // Number of discretization of different variables
        int nstep; // Number of discritization of the time step
        int itf;
        int n1;

        double ak10, t10; // t10 : conformal time till redshift 10
                          // ak10 : 500/Conftime_upto redshift 10
        int n10;          // Number of time grids after redshift 10

        int nriend;

        int itau; // Lensing Flag. itau = 1, lensing is required.

        int llo;
        double cllo, clhi;

        double dkn1, dkn2, dlnk0;
        int kkmin, kkmax; // index of minimum and maximum k in ak0[] array. Depends on user input. If someone don't want high k modes then kkmax will be at the middle. If low k modes are not required then kkmin will be in the middle.

        double *y, *yprime, *yt, *ytprime;

        // Work array for Dverk
        double **w, **wt;
        double *c;
        int ind;

        // Source Functions
        // Scalar
        double **d, **dp, **dkk, **phik;
        double **dadb, **dpadb, **dkkadb;
        double **diso, **dpiso, **dkkiso;

        double **dpr, **dppr, **dkpr;
        double **dpradb, **dppradb, **dkpradb;
        double **dpriso, **dppriso, **dkpriso;

        double phi; // Some temporary variable
            //Tensor
        double **dt, **dte, **dtb;
        double **dtpr, **dtepr, **dtbpr;

        // Source Function  // In New Grid for integration with bessel function
        // Scalar
        double *s2, *sp2, *sk2;
        double *s2adb, *sp2adb, *sk2adb;
        double *s2iso, *sp2iso, *sk2iso;

        // Tensor
        double *st2, *ste2, *stb2;

        //Brightness Fluctuation functions
        // For single field model
        double *dl2, *dpl2, *dkl2;
        //For twofield inflation
        double *dl2adb, *dpl2adb, *dkl2adb;
        double *dl2iso, *dpl2iso, *dkl2iso;

        // Brightness fluctuation functions
        double *dl3, *dpl3, *dkl3;
        double *dl3adb, *dpl3adb, *dkl3adb;
        double *dl3iso, *dpl3iso, *dkl3iso;

        double *dl, *dpl, *dkl;
        double *dladb, *dpladb, *dkladb;
        double *dliso, *dpliso, *dkliso;

        // Brightness fluctuations for tensor perturbations
        double *dtl, *dtel, *dtbl;
        double *dtl2, *dtel2, *dtbl2;
        double *dtl3, *dtel3, *dtbl3;

        // Power spectra for scalar perturbations before normalization
        double **cl, **cpl, **ccl, **ckkl, **ctkl;
        double **cladb, **cpladb, **ccladb;
        double **cliso, **cpliso, **ccliso;
        double **clcross, **cplcross, **cclcross;

        // Power spectra for tensor perturbations before normalization
        // Calculated for specific l
        double **ctl, **ctel, **ctbl, **ctcl;

        double **cltotal, **cpltotal, **ccltotal;
        double **clprtotal, **cplprtotal, **cclprtotal;

        double **clpr, **cplpr, **cclpr;
        double **ckklpr, **ctklpr;
        double **clpradb, **cplpradb, **cclpradb;
        double **clpriso, **cplpriso, **cclpriso;
        double **clprcross, **cplprcross, **cclprcross;

        double **ctlpr, **ctelpr, **ctblpr, **ctclpr;

        double d0lo = 1.0e40, d0hi = 1.0e40; // Spline Interpolation boundary condition

        int nw, nwt;

        int nk, nk2;

        double stpt = 50.0;

        double Dpass[4]; // For passing values - temp array
        double Dlpass[3];

        // Variables related to growth factor calculation
        double gftemp1, gftemp2, gftemp3, gftemp4, gftemp5, gftemp6, growthfactor;

        double *ak1, *dak1; // Wave number grid for integration of Source functions.

        int mxx[nstep0 + 1], m1, m2;

        int nstart1, nstop1, nstop1a;
        int nstart2, nstop2, nstop2a;

        FILE *fp[10];
        int printperturbations = 1;
        if (printperturbations == 1)
        {

            fp[1] = fopen("scalefactor.d", "w");
            fp[2] = fopen("ahdot.d", "w");
            fp[3] = fopen("eta.d", "w");
            fp[4] = fopen("deltac.d", "w");
            fp[5] = fopen("thetac.d", "w");
            fp[6] = fopen("deltab.d", "w");
            fp[7] = fopen("thetab.d", "w");
            fp[8] = fopen("deltagamma.d", "w");
            fp[9] = fopen("thetagamma.d", "w");
        }

        scalar Scalar;
        tensor Tensor;
        doubleinflation Doubleinflation;

        omegam = Variables.omegab + Variables.omegan + Variables.omegac; // Total matter density of the universe Baryonic + Neutrino + CDM

        // (P1) :
        // Calculating the final redshift at which the program will stop the calculations
        // arXiv:astro-ph/9603033
        // Seljak, Zaldarriaga (1996)  -  Page No 13. "Free streaming"

        if (fabs(Variables.omegab + Variables.omegac - 1.0) < .001 && Variables.zri == 0.0 && Variables.optdlss == 0.0 && Variables.itflag == 0 && Variables.h0 > 40.0) // The poin where the program stop calculation. If reionization or
            zst = 10.0;                                                                                                                                                 //	the tensor perturbations are requested then we should intigrate
        else                                                                                                                                                            // till the present era
            zst = 0.0;

        if (Variables.ict != 0)                                    // If transfer functions are requested then we should stop at the z
            zst = Numericx.min(zst, Variables.ztf[Variables.ntf]); // where the last transfer function is requested or after that

        // (P2) :
        if (Variables.itflag != 2)   // If scalar terms are present then we will take photon multipoles
            Variables.lmaxg = lmax0; // upto these many steps

        Variables.lmaxnr = lm2; //	We will take massless neutrino multipoles upto these many steps
        nstep = nstep0;         //	Number of discretization of the time step

        if (Variables.itflag != 0)    // If tensor perturbations are requested
            Variables.lmaxt = lmaxt0; // We will take photon multipoles upto these many steps
        else
            Variables.lmaxt = 0;

        // (P3) : Massive neutrino distribution

        /****************************************************************************************************************/
        /*		Initialize Massive Neutrino distribution										  										    	 */
        /*		"Cosmology"	- S.Weinberg																											 */
        /****************************************************************************************************************/

        if (Variables.annunr == 0 && Variables.omegan == 0.0) // If no massive neutrinos are present then set
        {                                                     // all the massive neutrino parameters to zero
            Variables.amnu = 0.0;
            Variables.nqmax = 0;
            Variables.lmaxnu = 0;
        }
        else // If neutrino present then
        {
            Variables.amnu = (Variables.omegan / Variables.annunr) * (Variables.grhom * 1000 / Variables.grhonr) * constan / (1.5 * zeta3); // Variables.amnu=m_nu*c**2/(k_B*T_nu0)
            // Probably its variables.amnu
            Variables.nqmax = nql;
            Variables.lmaxnu = lm3;

            double q, dq; // Temporary variables :)
            dq = 1.0;

            for (int i = 1; i <= Variables.nqmax; i++)
            {
                q = i * dq - 0.50;
                GlobalArray.dlfdlq[i] = -q / (1 + exp(-q)); // Calculate the d(log f)/d(log q)
            }
        }

        /*********************************************************************************************************************************/
        /*		Calculate number of equations                                                                                    */
        /*********************************************************************************************************************************/

        if (Variables.itflag != 2) // Scalers are requested
        {
            Variables.iq0 = 11 + 2 * Variables.lmaxg + Variables.lmaxnr; // Variables.iq0, Variables.iq1, Variables.iq2 will be used as an indicator
            Variables.iq1 = Variables.iq0 + Variables.nqmax;             // when we will store the massive neutrino number density,
            Variables.iq2 = Variables.iq1 + Variables.nqmax;             //	density and pressure, in function nu2()

            if (Variables.ndyn == 0 || Variables.ndyn == 3 || Variables.ndyn == 4)                                                  /********************************************************************/
            {                                                                                                                       /*  1) 3 				for storing a, a hdot, eta							  */
                Variables.nvar = 7 + 2 * (Variables.lmaxg + 1) + (Variables.lmaxnr + 1) + Variables.nqmax * (Variables.lmaxnu + 1); /*	 2) 2 				for storing CDM										  */
            }                                                                                                                       /*	 3) 2 				for storing Baryon									  */
                                                                                                                                    /*  4) 2 				for photon intensity and polarisation			  */
            else                                                                                                                    /*  5) 2*lmaxg 		        for storing photons higher multipoles in Legender Expansion							  */
            {                                                                                                                       /*	 6) lmaxnr       		for storing massless multipoles in Legender Expansion 					  */
                Variables.nvar = 9 + 2 * (Variables.lmaxg + 1) + (Variables.lmaxnr + 1) + Variables.nqmax * (Variables.lmaxnu + 1); /*	 7) 3*nqmax			for storing massive neutrinos						  */
            }                                                                                                                       /*	 8) nqmax*lmaxnu         	for Psi_nu for all Higher multipoles for massive neutrinos	  */
        }                                                                                                                           /*  9) 2 				for dark energy perturbation						  */
                                                                                                                                    /*	 																			  		  */
                                                                                                                                    /*	  This is how the above numbers came :)								  */
        else                                                                                                                        /********************************************************************/
        {
            Variables.nvar = 0; // If scalars are not required
        }

        if (Variables.itflag != 0) // Tensors are requested
        {
            Variables.nvart = 3 + 2 * (Variables.lmaxt + 1) + Variables.lmaxnr + 1; /*********************************************************************/
        }                                                                       /*	 1) 3					for a, h, h_dot											*/
        else                                                                    /*	 2) 2*(Variables.lmaxt+1)		for storing photon perturbations 									 	*/
        {                                                                       /*       3) (lmaxnr+1)                          for massless neutrinos									*/
            Variables.nvart = 0;                                                /*																							*/
        }                                                                       /*							*/

        /****************************************************************************************************************/
        /*		Initiallize Massive Neutrino distribution																				 		 */
        /****************************************************************************************************************/

        if ((Variables.annunr != 0) || (Variables.omegan != 0.0))
        {
            arel = 1.0 - 3.0 / Variables.amnu;
            Neutrino.initnul();
        }

        // (P4) : Calculating time of different epochs

        /****************************************************************************************************************/
        /*		Time at the present epoch																											 */
        /****************************************************************************************************************/
        Neutrino.initnul();

        Variables.tau0 = Numericx.rombint(Dtauda, 1.0e-8, 1.0, 1.0e-8); //tol); // This gives the conformal time in the present era
        Variables.epsw = 100.0 / Variables.tau0;                        // Used for calculating tight coupling in tensor perturbation

        /*****************************************************************************************************************/
        /*		Calculate the time of reionization																								  */
        /*****************************************************************************************************************/

        // (P5) : Calculating reionization

        if (Variables.optdlss > 0.0)
            Reion.reiopar(1.0); // It will calculate Variables.zri and Variables.zristp

        // (P6) :  k and t grid specification

        akmax = Variables.akmax0 / Variables.tau0; // Maximum and minimum k-values.
        akmin = 0.15 / Variables.tau0;

        if (Variables.itflag == 0)
            dlnk = 0.1; // dlnk is the logarithmic spacing used for low k.
        else            // if only scalars are requested then we large speacing can be taken
            dlnk = 0.05;

        if (akmax == 0)
            GlobalArray.dtaurec = 0.0; // This will only happen if there is some error in the user input for "ketamax"
        else
            GlobalArray.dtaurec = 4.0 / akmax;

        taumax = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / (zst + 1.0), Initval.tol); // taumax is the time where the program stops calculation

        dlntau0 = 0.0050;

        if (Variables.itflag != 0)
            dlntau0 = 0.0025; // if only tensor perturbations are required then take small step

        if (Variables.zri != 0.0)
        {
            Variables.taurist = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / (1 + Variables.zri), Initval.tol);     // Conformal time where the reionozation starts
            Variables.tauristp = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / (1 + Variables.zristp), Initval.tol); // Conformal tile where the reionization stops
        }
        else
        {
            Variables.taurist = Variables.tau0;
            Variables.tauristp = Variables.tau0;
        }

        if (Variables.ict == 0)
            taumin = 0.001 / akmax; // taumin if only Cl is requested
        if (Variables.ict == 1)
            taumin = 0.001 / Variables.akmaxt; // taumin if only TF is requested
        if (Variables.ict == 2)
            taumin = 0.001 * Numericx.min(1.0 / akmax, 1.0 / Variables.akmaxt); //	taumin if both are needed

        taumin = Numericx.min(taumin, 0.1); // Start in the radiation dominated era

        if (Variables.amnu != 0.0)
            taumin = Numericx.min(taumin, arel / Variables.adotrad); // Start when massive neutrinos were strongly relativistic

        /*************************************************************************************************************/
        /* 			Initialize baryon temperature and ionization fractions vs. time.											 */
        /*************************************************************************************************************/
        double tautemp[2], dlsso;
        int nsteptemp[1];
        nsteptemp[0] = nstep;                                                       //
        n1 = Recombination.finithermo(taumin, taumax, dlntau0, nsteptemp, tautemp); //output - > taurend, n1,
        nstep = nsteptemp[0];
        taurend = tautemp[0];

        // (P7) : k splitting

        dlsso = Variables.tau0 - Variables.taurmax;

        shor = Numericx.rombint(Dsoundda, 1.0e-8, GlobalArray.armax, Initval.tol); // Sound Horizon at LSS

        Variables.aksplit = (1.5 / shor) * Variables.aksplit; // K splitting
        // Compute cut wavevector

        if ((Variables.kcutflag == 1) || (Variables.kcutflag == -1)) // Output dlss to shift spectra.
        {
            printf("### First line has dlss for this model and kcutflag"); // Write dlss in the output file.
            printf("### Then normal l ClT ClE ClC output");
            printf("%e %d", dlsso, Variables.kcutflag);
        }

        // (P8) : Calculating the points where the transfer functions will be stored

        if (Variables.itflag != 2) // Only tensors are requested
        {
            d = matrix(nk0 + 1, nstep0 + 1); // For the source functions
            dp = matrix(nk0 + 1, nstep0 + 1);
            dkk = matrix(nk0 + 1, nstep0 + 1);
            phik = matrix(nk0 + 1, ntau + 1);

            dpr = matrix(nk0 + 1, nstep0 + 1); // For interpolating the Source Functions in the integration grid
            dppr = matrix(nk0 + 1, nstep0 + 1);
            dkpr = matrix(nk0 + 1, nstep0 + 1);

            if (Variables.twofield == 1) // Two fields inflation is requested
            {
                dadb = matrix(nk0 + 1, nstep0 + 1);
                dpadb = matrix(nk0 + 1, nstep0 + 1);
                dkkadb = matrix(nk0 + 1, nstep0 + 1);

                diso = matrix(nk0 + 1, nstep0 + 1);
                dpiso = matrix(nk0 + 1, nstep0 + 1);
                dkkiso = matrix(nk0 + 1, nstep0 + 1);

                dpradb = matrix(nk0 + 1, nstep0 + 1);
                dppradb = matrix(nk0 + 1, nstep0 + 1);
                dkpradb = matrix(nk0 + 1, nstep0 + 1);

                dpriso = matrix(nk0 + 1, nstep0 + 1);
                dppriso = matrix(nk0 + 1, nstep0 + 1);
                dkpriso = matrix(nk0 + 1, nstep0 + 1);
            }
        }

        if (Variables.itflag != 0) //If tensors are required
        {
            dt = matrix(nk0 + 1, nstep0 + 1);
            dte = matrix(nk0 + 1, nstep0 + 1);
            dtb = matrix(nk0 + 1, nstep0 + 1);

            dtpr = matrix(nk0 + 1, nstep0 + 1);
            dtepr = matrix(nk0 + 1, nstep0 + 1);
            dtbpr = matrix(nk0 + 1, nstep0 + 1);
        }

        tautf0 = (double *)malloc((ntfmax + 1) * sizeof(double));

        if (Variables.ict != 0) // Transfer functions are required
        {
            a0 = 1.0e-8;
            tautf[0] = 0.0;
            aprevious = a0;

            double atf;
            for (itf = 1; itf <= Variables.ntf; itf++) // Calculating the times for the outputs of the transfer functions.
            {
                atf = 1.0 / (Variables.ztf[itf] + 1.0);
                tautf[itf] = tautf[itf - 1] + Numericx.rombint(Dtauda, aprevious, atf, Initval.tol);
                tautf[itf] = Numericx.min(tautf[itf], GlobalArray.atau0[nstep]);
                indtau[itf] = itf;
                aprevious = atf;
            }

            if (Variables.lensflag != 0)
            {
                for (itf = 1; itf <= Variables.ntf; itf++)
                    tautf0[itf] = tautf[itf];

                for (int i = 1; i <= ntau + Variables.ntf; i++)
                    indtau[i] = 0;

                double ar = 1.0 / 1090.0;
                Variables.taur = Numericx.rombint(Dtauda, a0, ar, Initval.tol);
                itf = 1;

                double taunew = Variables.taur;

                for (int i = 1; i <= ntau; i++)
                {
                    if (taunew <= tautf0[itf] || itf > Variables.ntf)
                        tautf[i + itf - 1] = taunew;

                    else
                    {
                        tautf[i + itf - 1] = tautf0[itf];
                        indtau[i + itf - 1] = itf;
                        itf = itf + 1;
                        continue;
                    }
                    taunew = i * (Variables.tau0 - Variables.taur) / (ntau - 1) + Variables.taur;
                }

                if (itf <= Variables.ntf)
                {
                    for (int i = itf; i <= Variables.ntf; itf++)
                    {
                        indtau[ntau + i] = i;
                        tautf[ntau + i] = tautf0[i];
                    }
                }
                Variables.ntf = Variables.ntf + ntau;
            }
        }
        tautf[Variables.ntf - 1] = tautf[Variables.ntf - 1] - Initval.tol;
        free(tautf0);

        // Integration will be carried out after z=10 only for low
        // k, k < k10. If there is reionization the boundary will not be z=10 but Variables.tauristp

        // No need to calculate for the large wavelengths after z=10

        if (zst != 10.0)
        {
            t10 = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / 11.0, Initval.tol);

            ak10 = 500.0 / Variables.tau0;

            n10 = n1 + int(log(t10 / taurend) / dlntau0);

            if (Variables.zri != 0.0)
            {
                nriend = GlobalArray.j2ri1 + GlobalArray.nri + GlobalArray.nri0;
                t10 = Numericx.max(t10, GlobalArray.atau0[nriend]);
                n10 = nriend + int(log(t10 / GlobalArray.atau0[nriend]) / dlntau0);
            }
        }
        else
        {
            t10 = Variables.tau0;
            ak10 = akmax;
            n10 = nstep + 1;
        }

        // (P9): 
        /************************************************/
        /*          If we want to calculate Cl          */
        /************************************************/
        ak0 = (double *)malloc((nk0 + 1) * sizeof(double)); // Allocate the array for k grid points

        if (Variables.ict != 1)
        {
            /*************************************************************************************************************************/
            /* set k values for which the sources for the anisotropy and polarization will be calculated. 									 */
            /* For low values of k we use a logarithmic spacing.																							 */
            /*************************************************************************************************************************/

            if (Variables.zri != 0.0)
                dlnk0 = 2.0 * dlnk;
            else
                dlnk0 = 5.0 * dlnk; // if no reionization then take larger steps in the calculation

            dkn1 = 0.8 / Variables.taurmax;
            dkn2 = 1.5 / Variables.taurmax;

            int nk1 = int(log(dkn1 / (akmin * dlnk0)) / dlnk0) + 1;

            double akchange = 5.0 * 3.14159265 / shor;

            if (akmax > akchange)
            {
                nk2 = int((akchange - akmin * exp((nk1 - 1) * dlnk0)) / dkn1) + nk1 + 1;
                nk = int((akmax - akchange) / dkn2) + nk2 + 1;
            }
            else
            {
                nk = int((akmax - akmin * exp((nk1 - 1) * dlnk0)) / dkn1) + nk1 + 1;
                nk2 = nk + 1;
            }

            if (nk > nk0)
            {
                printf("Sorry, the arrays were dimensioned for a maximum of");
                printf("%d, k modes. The model you requested needs, %d", nk0, nk);
                printf("Please make the arrays bigger by making ");
                printf("nk0 bigger where it appears");
                exit(1);
            }

            for (int ik = 1; ik <= nk; ik++)
            {
                if (ik <= nk1)
                    ak0[ik] = akmin * exp((ik - 1) * dlnk0);
                else
                {
                    if (ik > nk2)
                        ak0[ik] = ak0[nk2] + (ik - nk2) * dkn2;
                    else
                        ak0[ik] = ak0[nk1] + (ik - nk1) * dkn1;
                }
            }

            kkmin = 1;
            kkmax = nk;

            /************************************************************************************************************************/
            /*    K split																																				*/
            /************************************************************************************************************************/
            if ((Variables.kcutflag) == 1)
            {
                for (int k = 1; k <= nk; k++)
                {
                    if (ak0[k] < (1.3 * Variables.aksplit))
                        kkmax = k;
                }

                kkmax = Numericx.minint(nk, kkmax + 4);
            }

            if ((Variables.kcutflag) == -1)
            {
                for (int k = 1; k <= nk; k++)
                {
                    if (ak0[k] < (0.3 * Variables.aksplit))
                        kkmin = k;
                }
                kkmin = Numericx.maxint(1, kkmin - 4);
            }

            /*********************************************************************************************************************************/
            /*  Main calculation. Loop over wavenumbers.																													*/
            /*  For calculation of Cl we have to calculate the source terms in the different time and different wavenumber							*/
            /*********************************************************************************************************************************/

            if (Variables.itflag != 2) // If scalars are requested
            {
                y = (double *)malloc((nvar0 + 1) * sizeof(double));
                yprime = (double *)malloc((nvar0 + 1) * sizeof(double));
            }
            if (Variables.itflag != 0) // If tensors are requested
            {
                yt = (double *)malloc((nvar0t + 1) * sizeof(double));
                ytprime = (double *)malloc((nvar0t + 1) * sizeof(double));
            }
        L002:
            //            for (int iii = 1; iii <= Variables.nvar; iii++)
            //                yprime[iii] = 0.0;

            c = (double *)malloc(25 * sizeof(double));

            w = (double **)malloc(10 * sizeof(double));
            for (int i = 0; i < 10; i++)
                w[i] = (double *)malloc((nvar0 + 1) * sizeof(double));

            wt = (double **)malloc(10 * sizeof(double));
            for (int i = 0; i < 10; i++)
                wt[i] = (double *)malloc((nvar0t + 1) * sizeof(double));

            for (int ik = kkmin; ik <= kkmax; ik++) // Minimum and maximum value of the webnumber
            {
                Variables.ak = ak0[ik]; // Discretized values of k's.
                // Begin when wave is far outside horizon.
                taustart = 0.001 / Variables.ak; // Conformal time (in Mpc) in the radiation era.
                // Start calculations early in the radiation era.
                taustart = Numericx.min(taustart, 0.1); // We are not including neutrinos as sources of the tensor modes.

                taustart0 = taustart; // Starting time should be saved as required for further calculations

                if (Variables.amnu != 0.0) // Start when massive neutrinos are strongly relativistic.
                {
                    arel = 1 - 3 / Variables.amnu;
                    taustart = Numericx.min(taustart0, arel / Variables.adotrad); // In radiation era a ~ t
                }
                Variables.ak2 = Variables.ak * Variables.ak; // k^2 = k*k;
                if (Variables.itflag != 2)
                {
                    // Variables.itflag = 2 -> Only tensors are requested
                    Scalar.finitial(y, taustart); // As Cl is required so initilize the terms for calculating the source

                    if (printperturbations == 1)
                    {
                        fprintf(fp[1], "%e  ", y[1]);
                        fprintf(fp[2], "%e  ", y[2]);
                        fprintf(fp[3], "%e  ", y[3]);
                        fprintf(fp[4], "%e  ", y[4]);
                        fprintf(fp[5], "%e  ", y[5]);
                        fprintf(fp[6], "%e  ", y[6]);
                        fprintf(fp[7], "%e  ", y[7]);
                        fprintf(fp[8], "%e  ", y[8]);
                        fprintf(fp[9], "%e  ", y[9]);

                        //                        printf("%e  \n",y[1]);
                        //                        printf("%e  \n",y[2]);
                        //                        printf("%e  \n",y[3]);
                        //                        printf("%e  \n",y[4]);
                        //                        printf("%e  \n",y[5]);
                        //                        printf("%e  \n",y[6]);
                        //                        printf("%e  \n",y[7]);
                        //                        printf("%e  \n",y[8]);
                        //                        printf("%e  \n",y[9]);

                        //                       for(int iuiu; iuiu<30;iuiu++)
                        //                       {
                        //                           printf("%d  %e \n",iuiu,y[iuiu]);
                        //                      }
                    }

                    if (Variables.ndyn == 1 || Variables.ndyn == 2) // ndyn =1,2 -> Dark energy is perturbated
                    {
                        y[Variables.nvar - 1] = y[nvar0 - 1];
                        y[Variables.nvar] = y[nvar0];
                    }
                    tau = taustart; // d[] -> the sources for the anisotropy,
                    // dp[] -> sources for the polarization,   t stands for tensor.
                    if (Variables.ict != 0)
                        itf = 1; //	Transfer function is required  [transfer function loov var  1 <= itf <= nft ]
                    if (Variables.lensflag != 0)
                        itau = 1; // Lensing is required

                    double tol1 = Numericx.min(Initval.tol, (1.0e-17) * (Variables.ak2 / 1.0e-10)); // Tollerance value. The instability can ocure at low k. So increase the tolerance at low k

                    if (Variables.ict == 0) // No transfer functions are requested
                        taulasttf = 0.0;
                    else
                        taulasttf = tautf[Variables.ntf]; //  It stores the time time where the last transfer function is requested

                    /******************************************************************************************************************************/
                    /*  Integrator Communication Variables																														*/
                    /******************************************************************************************************************************/

                    ind = 2; // This is specified for the "dverk" intigrator
                    for (int j = 1; j <= 24; j++)
                        c[j] = 0.0; // This is specified for the "dverk" intigrator

                    c[3] = 1.0e-8; // This is specified for the "dverk" intigrator

                    nw = nvar0;


                    /******************************************************************************************************************************/
                    /*  Loop over timesteps.																																		*/
                    /*	 This part only calculate the Scaler Source Terms and the Transfer Function																*/
                    /******************************************************************************************************************************/
                    for (int j = 2; j <= nstep; j++)
                    {
                        tauend = GlobalArray.atau0[j];
                        if (Variables.ak > ak10 && tauend > t10 && (Variables.ict == 0 || tau > tautf[Variables.ntf])) // ak10 & t10 srtores maximum values of k and t to be intigrated
                        {                                                                                              // So if k or t crosses this limit then just store 0 in the source terms
                            d[ik][j] = 0.0;
                            dp[ik][j] = 0.0;
                        }
                        else
                        {
                            {
                                tau = Numericx.dverk(Variables.nvar, Fderivs, tau, y, tauend, tol1, ind, c, nw, w); // Intigrator
                                ind = 3;
                                //printf("%e \n",tau);
                            }
                            if (printperturbations == 1 && j < 2550)
                            {
                                fprintf(fp[1], "%e  ", y[1]);
                                fprintf(fp[2], "%e  ", y[2]);
                                fprintf(fp[3], "%e  ", y[3]);
                                fprintf(fp[4], "%e  ", y[4]);
                                fprintf(fp[5], "%e  ", y[5]);
                                fprintf(fp[6], "%e  ", y[6]);
                                fprintf(fp[7], "%e  ", y[7]);
                                fprintf(fp[8], "%e  ", y[8]);
                                fprintf(fp[9], "%e  ", y[9]);

                                //printf("%e  ",y[1]);
                                //printf("%e  ",y[2]);
                                //printf("%e  ",y[3]);
                                //printf("%e  ",y[4]);
                                //printf("%e  ",y[5]);
                                //printf("%e  ",y[6]);
                                //printf("%e  ",y[7]);
                                //printf("%e  ",y[8]);
                                //printf("%e  ",y[9]);
                            }
                            Scalar.fderivs(Variables.nvar, tau, y, yprime);
                            Scalar.foutput(Variables.nvar, y, yprime, j, Variables.tau0, tau, Dpass);

                            d[ik][j] = Dpass[0];
                            dp[ik][j] = Dpass[1];
                            dkk[ik][j] = Dpass[2];
                            phi = Dpass[3];

                            /************************************************************************************************************************/
                            /*  Calculate the transfer function.																												*/
                            /************************************************************************************************************************/
                            do //for(;(j<nstep && itf<=ntf && GlobalArray.atau0[j+1]>tautf[itf]);)
                            {
                                if (Variables.ict != 0 && itf <= Variables.ntf) // IF NUMBER J							// itc =0 -> Only Power spectrum, ntf -> number of transfer function required
                                {

                                    if (j < nstep && tauend < tautf[itf] && GlobalArray.atau0[j + 1] > tautf[itf])
                                    {
                                        tau = Numericx.dverk(Variables.nvar, Fderivs, tau, y, tautf[itf], tol1, ind, c, nw, w); // Intigrator
                                        ind = 3;
                                    }
                                    // output transfer functions for this k-value.
                                    if (fabs(tau - tautf[itf]) < 1.0e-5) // IF NUMBER K
                                    {
                                        if (indtau[itf] != 0)
                                            Scalar.outtransf(Variables.nvar, y, 0.0, indtau[itf]); // It will calculate and store the transfer function
                                        else
                                        {
                                            if (Variables.lensflag != 0)
                                            {
                                                Scalar.fderivs(Variables.nvar, tau, y, yprime);
                                                Scalar.foutput(Variables.nvar, y, yprime, j, Variables.tau0, tau, Dpass); // For lensing we only need the lensing parameter phi. so store the sources in temporary variable
                                                phi = Dpass[3];
                                                phik[ik][itau] = phi;
                                                itau = itau + 1;
                                            } // IF END
                                        }     // ELSE END

                                        itf = itf + 1;
                                    }                                                                                                           // IF NUMBER K END
                                }                                                                                                               // IF NUMBER J END
                            } while ((j < nstep) && (itf <= Variables.ntf) && (GlobalArray.atau0[j + 1] > tautf[itf]) && (Variables.ict != 0)); // FOR ENDS HERE
                        }                                                                                                                       // else end here
                        //					fprintf(fdkk,"%e\t",dkk[ik][j]);
                    } // FOR OVER J WILL END HERE

                    //fprintf(fdkk,"\n");

                    if (printperturbations == 1)
                    {
                        fprintf(fp[1], "\n");
                        fprintf(fp[2], "\n");
                        fprintf(fp[3], "\n");
                        fprintf(fp[4], "\n");
                        fprintf(fp[5], "\n");
                        fprintf(fp[6], "\n");
                        fprintf(fp[7], "\n");
                        fprintf(fp[8], "\n");
                        fprintf(fp[9], "\n");
                    }

                    d[ik][nstep] = 0.0;
                    dp[ik][nstep] = 0.0;
                } // IF OVER Variables.itflag ENDS HERE

                /******************************************************************************************************************************/
                /*  Loop over timesteps.																																		*/
                /*	 Tensors will only be calculated if k*tau< stsp.																									*/
                /******************************************************************************************************************************/

                stpt = 50.0;

                if (Variables.itflag != 0)
                {
                    Tensor.finitialt(taustart0, yt);
                    tau = taustart0;

                    /******************************************************************************************************************************/
                    /*  Integrator Communication Variables																														*/
                    /******************************************************************************************************************************/

                    ind = 1; // This is specified for the "dverk" integrator
                    nwt = Variables.nvart + 5;
                    /***************************************************************************************************************************/
                    /*  Loop over timesteps.																																	*/
                    /*	 dt temp, dte electric field, dtb magnetic field		 																						*/
                    /***************************************************************************************************************************/
                    Tensor.fderivst(Variables.nvart, tau, yt, ytprime);
                    for (int j = 2; j <= nstep; j++)
                    {
                        tauend = GlobalArray.atau0[j];
                        if ((Variables.ak * tauend) > stpt)
                        {
                            dt[ik][j] = 0.0;
                            dte[ik][j] = 0.0;
                            dtb[ik][j] = 0.0;
                        }
                        else
                        {
                            double tol1 = Initval.tol;
                            {
                                tau = Numericx.dverk(Variables.nvart, Fderivst, tau, yt, tauend, tol1, ind, c, nwt, wt); // Intigrator
                                ind = 3;
                            }

                            //for(int abc=0;abc<10;abc++)
                            //    printf("%e  ",yt[abc]);
                            //printf("\n");

                            Tensor.fderivst(Variables.nvart, tau, yt, ytprime);                            // Tensor derivatives
                            Tensor.foutputt(Variables.nvart, yt, ytprime, j, Variables.tau0, tau, Dlpass); // Tensor source perturbation
                            dt[ik][j] = Dlpass[0];
                            dte[ik][j] = Dlpass[1];
                            dtb[ik][j] = Dlpass[2];
                        }
                    }
                    dt[ik][nstep] = 0.0;
                    dte[ik][nstep] = 0.0;
                    dtb[ik][nstep] = 0.0;
                }

            } // FOR LOOP OVWE WAVE NUMBER ENDS HERE

            /*****************************************************************************************************************************/
            //	This part is added for the two field inflation models 																						  //
            /*****************************************************************************************************************************/

            free(c);
            free(w);
            free(wt);

            if (Variables.itflag != 0) // If tensors are requested
            {
                free(yt);
                free(ytprime);
            }
            if (Variables.itflag != 2) // If scalars are requested
            {
                free(y);
                free(yprime);
            }

            if (Variables.twofield == 1)
            {
                for (int storei = 0; storei < nk0 + 1; storei++)
                    for (int storej = 0; storej < nstep0 + 1; storej++)
                    {
                        if (Variables.initfl == 1)
                        {
                            dadb[storei][storej] = d[storei][storej];
                            dpadb[storei][storej] = dp[storei][storej];
                            dkkadb[storei][storej] = dkk[storei][storej];
                        }
                        else if (Variables.initfl == 2)
                        {
                            diso[storei][storej] = d[storei][storej];
                            dpiso[storei][storej] = dp[storei][storej];
                            dkkiso[storei][storej] = dkk[storei][storej];
                        }
                        else
                        {
                            printf("\nERROR: In two field inflation only CDM isocurvature models are expected.");
                            exit(1);
                        }
                    }
                if (Variables.initfl == 1)
                {
                    Variables.initfl = 2;
                    goto L002; //	Redo the entire thing for CDM Isocurvature modes
                }
            }

            /*****************************************************************************************************************************/
        } // IF Variables.ict !=1 ENDS HERE

        if (Variables.ict == 1)
            nk = 0; // It is set to zero because it will be required later if only transfer functions are requested only.
        // (P10) : If Transfer functions are requested
        printf("\nIntegration over the perturbation variables for calculating Cl are finished\n\n");
        /******************************************************************************************************************************/
        /*  If the transfer functions are requested then Variables.ict!=0. So the transfer functions are calculated here	  							*/
        /******************************************************************************************************************************/
        int j;
        
        if (Variables.ict != 0) // If transfer functions are requested
        {
            int nkt;

            int closestknum;

            double akdone, akgf; // k-growth factor
            double phi1, phihk;
            if (Variables.ict == 2)
            {
                if (Variables.akmaxt > akmax)
                {
                    nkt = nk + int((log(Variables.akmaxt) - log(akmax)) * Variables.nlnkt) + 1;
                    akdone = ak0[nk];
                }
                else
                    nkt = nk;
            }
            else
            {
                nkt = int((log(Variables.akmaxt) - log(akmin)) * Variables.nlnkt) + 1;
                akdone = akmin;
            }

            akgf = 1.0;

            if (akdone > akgf)
                akgf = akdone * exp(1.5 / Variables.nlnkt);

            /*****************************************************************************************************************************/
            /*	 Loop over wavenumbers.																																	  */
            /*****************************************************************************************************************************/

            double exptemp;

            // Allocate dverk variables
            // ---------------------------------------------------
            c = (double *)malloc(25 * sizeof(double));

            w = (double **)malloc(10 * sizeof(double));
            for (int i = 0; i < 10; i++)
                w[i] = (double *)malloc((nvar0 + 1) * sizeof(double));

            y = (double *)malloc((nvar0 + 1) * sizeof(double));
            yprime = (double *)malloc((nvar0 + 1) * sizeof(double));
            for (int ik = nk + 1; ik <= nkt; ik++)
            {
                exptemp = (double)(ik - nk) / Variables.nlnkt;

                Variables.ak = akdone * exp(exptemp); // Calculate k, steps are exponentially seperated.
                //printf("\n exp = %e",exp(exptemp));
                // |...................|...........|
                //				      done(nk) 		(nkt)
                if (Variables.lensflag != 0)
                    itau = 1;

                /******************************************************************************************************************************/
                /*  Integrator Communication Variables																														*/
                /******************************************************************************************************************************/
                ind = 2;

                for (int jt = 1; j <= 24; j++)
                    c[jt] = 0.0;

                c[3] = 1.0e-10;

                /******************************************************************************************************************************/
                /*	 Time loop																																						*/
                /*  Begin when wave is far outside horizon.																												*/
                /*	 Conformal time (in Mpc) in the radiation era, for photons plus 3 species of relativistic neutrinos.								*/
                /******************************************************************************************************************************/
                taustart = 0.001 / Variables.ak;

                taustart = Numericx.min(taustart, 0.1); // Start early in the radiation era.

                if (Variables.amnu != 0.0)
                { // Start when massive neutrinos are strongly relativistic.
                    arel = 1.0 - 3 / Variables.amnu;
                    taustart = Numericx.min(taustart, arel / Variables.adotrad);
                }

                Variables.ak2 = Variables.ak * Variables.ak; // ak2 = k^2
                Scalar.finitial(y, taustart);                // initiallize source terms in y at the time taustart

                tau = taustart;

                j = nk;

                if (Variables.ak < akgf)
                {
                    for (int itf = 1; itf <= Variables.ntf; itf++)
                    {
                        if (fabs(tautf[itf] - tau) > 1.0e-5)
                        {
                            tau = Numericx.dverk(Variables.nvar, Fderivs, tau, y, tautf[itf], Initval.tol, ind, c, nw, w); // Intigrator
                            ind = 3;
                        }
                        else
                            tau = tautf[itf];

                        if (indtau[itf] != 0) // If transferfunction requied at this point then calculate that.
                            Scalar.outtransf(Variables.nvar, y, 0.0, indtau[itf]);
                        else
                        {
                            if (Variables.lensflag != 0) // If lensing required then store phi
                            {
                                if (itf == 1)
                                {
                                    Scalar.fderivs(Variables.nvar, tau, y, yprime);
                                    Scalar.foutput(Variables.nvar, y, yprime, j, Variables.tau0, tau, Dpass);
                                    phi1 = Dpass[3];
                                }

                                Scalar.fderivs(Variables.nvar, tau, y, yprime);
                                Scalar.foutput(Variables.nvar, y, yprime, j, Variables.tau0, tau, Dpass);
                                phi = Dpass[3];

                                phik[ik][itau] = phi;
                                printf("%e   ", phi);
                                itau = itau + 1;
                            }
                        }
                    }

                    closestknum = GlobalArray.istore[1] - 1;
                }
                else //  If k > final k specified then just store the powers for TF
                {

                    {
                        tau = Numericx.dverk(Variables.nvar, Fderivs, tau, y, tauend, Initval.tol, ind, c, nw, w); // Intigrator
                        ind = 3;
                    }

                    if (indtau[itf] != 0) //  If transfer function is requested at this point then call outtrans()
                        Scalar.outtransf(Variables.nvar, y, 0.0, indtau[itf]);

                    for (int itf = 2; itf <= Variables.ntf; itf++)
                    {
                        if (indtau[itf] != 0)
                        {
                            gftemp1 = Variables.omegac * GlobalArray.power[6 * closestknum + 2][indtau[itf]]; //
                            gftemp2 = Variables.omegab * GlobalArray.power[6 * closestknum + 3][indtau[itf]]; //
                            gftemp3 = Variables.omegan * GlobalArray.power[6 * closestknum + 6][indtau[itf]]; //
                            gftemp4 = Variables.omegac * GlobalArray.power[6 * closestknum + 2][1];           //
                            gftemp5 = Variables.omegab * GlobalArray.power[6 * closestknum + 3][1];           //	  Growth function calculation.
                            gftemp6 = Variables.omegan * GlobalArray.power[6 * closestknum + 6][1];           //	  Dodelson P.183
                            //
                            growthfactor = (gftemp1 + gftemp2 + gftemp3) / (gftemp4 + gftemp5 + gftemp6); //

                            GlobalArray.power[6 * GlobalArray.istore[indtau[itf] + 1]][itf] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 1][1]; //  Power
                            //  All followings are transfer functions
                            GlobalArray.power[6 * GlobalArray.istore[itf] + 2][indtau[itf]] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 2][1] * growthfactor;         //
                            GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 3][indtau[itf]] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 3][1] * growthfactor; //
                            GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 4][indtau[itf]] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 4][1] * growthfactor; //
                            GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 5][indtau[itf]] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 5][1] * growthfactor; //
                            GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 6][indtau[itf]] = GlobalArray.power[6 * GlobalArray.istore[indtau[itf]] + 6][1] * growthfactor; //
                            GlobalArray.istore[indtau[itf]] = GlobalArray.istore[indtau[itf]] + 1;
                        }
                        else
                        {
                            if (Variables.lensflag != 0)
                            {
                                Scalar.fderivs(Variables.nvar, tau, y, yprime);
                                Scalar.foutput(Variables.nvar, y, yprime, j, Variables.tau0, tau, Dpass);
                                phihk = Dpass[3];
                                phik[ik][itau] = phihk * phi / phi1;
                                itau = itau + 1;
                            } // END IF
                        }     // END ELSE
                    }         // END FOR
                }             // END ELSE
            }                 // END FOR WAVE NUMBER LOOP
            
            free(y);
            free(yprime);

            free(c);
            free(w);
        } // IF OF Variables.ict!=0 ENDS HERE

        if (Variables.lensflag != 0)
            Variables.ntf = Variables.ntf - ntau;

        if (Variables.ict != 0)
            output_power(Variables.ntf, Variables.amnu);

        // (P13) : Calculating the Cl
        /*********************************************************************************************************************************/
        /*  If CMB calculations are requested, calculate the Cl by integrating the sources over time and over k. 								*/
        /*********************************************************************************************************************************/
        if (Variables.ict != 1)
        {

            /******************************************************************************************************************************/
            /*  If Scaler Cl's are requested then initiallize the Cls and calculate the derivative of the source terms.							*/
            /******************************************************************************************************************************/

            //for (int ik = kkmin; ik <= kkmax; ik++)
            //    for (int j = 2; j <= nstep; j++)
            //    {
            // printf("%e\t", dpr[ik][j]);
            //        }

            if (Variables.itflag != 2)
            {
                if (Variables.twofield == 0)
                {
                    Numericx.spline2D(ak0, d, nk, nstep, d0lo, d0hi, dpr);    //
                    Numericx.spline2D(ak0, dp, nk, nstep, d0lo, d0hi, dppr);  // Calculate the double time derivative of the
                    Numericx.spline2D(ak0, dkk, nk, nstep, d0lo, d0hi, dkpr); // scaler source terms

                    /*
                    FILE *fpDp;
                    fpDp = fopen("Dp.txt","w");
                    FILE *fpDpr;
                    fpDpr = fopen("Dpr.txt","w");
                    for (int jtj = 1; jtj < nstep0 + 1; jtj ++){
                         for (int iti = 1; iti < nk0 + 1; iti++){
                            fprintf(fpDp,"\t%e", d[iti][jtj]);
                            fprintf(fpDpr,"\t%e",dpr[iti][jtj]);
                         }
                            fprintf(fpDp,"\n");
                            fprintf(fpDpr,"\n");                         
                    }
                    fclose(fpDp);
                    fclose(fpDpr);
                    exit(1);*/
                }
                /***************************************************************************************************************************************/
                // This part is added for two field inflation 																														//
                /***************************************************************************************************************************************/
                else if (Variables.twofield == 1) // For two field model
                {
                    Numericx.spline2D(ak0, dadb, nk, nstep, d0lo, d0hi, dpradb);    //
                    Numericx.spline2D(ak0, dpadb, nk, nstep, d0lo, d0hi, dppradb);  // Calculate the time derivative of the
                    Numericx.spline2D(ak0, dkkadb, nk, nstep, d0lo, d0hi, dkpradb); // scaler source terms

                    Numericx.spline2D(ak0, diso, nk, nstep, d0lo, d0hi, dpriso);    //
                    Numericx.spline2D(ak0, dpiso, nk, nstep, d0lo, d0hi, dppriso);  // Calculate the time derivative of the
                    Numericx.spline2D(ak0, dkkiso, nk, nstep, d0lo, d0hi, dkpriso); // scaler source terms
                }
            }
            //            for (int ik = kkmin; ik <= kkmax; ik++)
            //                for (int j = 2; j <= nstep; j++)
            //                {
            //                    printf("%d %d %e %e %e %e %e %e\n", ik, j, d[ik][j], dp[ik][j], dkk[ik][j], dpr[ik][j], dppr[ik][j], dkpr[ik][j]);
            //                    }
            /******************************************************************************************************************************/
            /*  If Tensor Cl's are requested then initiallize the Cls and calculate the derivative of the source terms.							*/
            /******************************************************************************************************************************/

            if (Variables.itflag != 0)
            {
                Numericx.spline2D(ak0, dt, nk, nstep, d0lo, d0hi, dtpr);   //
                Numericx.spline2D(ak0, dte, nk, nstep, d0lo, d0hi, dtepr); // Calculate the time derivative of the
                Numericx.spline2D(ak0, dtb, nk, nstep, d0lo, d0hi, dtbpr); // tensor source terms
            }

            // (P14) : Define the wavenumber and the difference grid for integration
            /******************************************************************************************************************************/
            /*  Fixing the wavenumbers for integration.																												*/
            /******************************************************************************************************************************/

            double dk = 2.5 / Variables.tau0;
            double dk0 = 1.5 / Variables.tau0;
            int no = 700; // The dimension of K arrays rellated to no =700
            double dlnk1 = 0.07;

            int nko, no1;

            no1 = (int)(log(10.0 * dk0 / akmin) / dlnk1) + 1;

            if (akmax > (no * dk0))
                nko = (int)((akmax - no * dk0) / dk) + no;
            else
            {
                no = (int)((akmax - 10.0 * dk0) / dk0) + no1;
                nko = no;
            }

            if (nko > nkmax) // Dimension of K array is nkmax
            {
                printf("\nSorry, the arrays were dimensioned for a maximum of %d k modes", nkmax);
                printf("\nThe model you requested needs %d", nko);
                printf("\nPlease make the arrays bigger by making nkmax bigger where it appears");
                exit(1);
            }

            ak1 = (double *)malloc((nkmax + 1) * sizeof(double));
            dak1 = (double *)malloc((nkmax + 1) * sizeof(double));

            for (int k = 1; k <= nko; k++)                             //   |...................|...............|..........................|
            {                                                          // 	   <for short wev> (no1) <mid wave> (no)  <long wave length>
                if (k <= no)                                           //        log space            linear         linear large space
                {                                                      //
                    if (k <= no1)                                      //    Wave number grid points are stored in ak[]
                    {                                                  //    Wave number spacing are stoder in dak1[]
                        ak1[k] = 10.0 * dk0 * exp(-(no1 - k) * dlnk1); //
                        dak1[k] = ak1[k] * dlnk1;                      //
                    }
                    else
                    {
                        ak1[k] = ak1[no1] + (k - no1) * dk0;
                        dak1[k] = dk0;
                    }
                }
                else
                {
                    ak1[k] = ak1[no] + (k - no) * dk;
                    dak1[k] = dk;
                }
            }

            dak1[1] = 0.5 * dak1[1];
            dak1[no1] = 0.5 * (dak1[no1] + dk0);
            dak1[nko] = 0.5 * dak1[nko];

            /*****************************************************************************************************************************/
            /*   K split																																					  */
            /*****************************************************************************************************************************/

            // |Kmin..............................|Kmax

            kkmin = 1;   // Minimum Value of K (Wavenumber)
            kkmax = nko; // Maximum value of K (Wavenumber)

            if (Variables.kcutflag == 1)       // |Kmin............|Kmax.................|
            {                                  // |1...............|1.3..................|
                for (int k = 1; k <= nko; k++) // K < K*
                    if (ak1[k] < 1.3 * Variables.aksplit)
                        kkmax = k + 1;
            }

            if (Variables.kcutflag == -1)      // |............|Kmin.................|Kmax
            {                                  // |............|0.3..................|nko
                for (int k = 1; k <= nko; k++) // K > K*
                    if (ak1[k] <= 0.3 * Variables.aksplit)
                        kkmin = k;
            }

            /*****************************************************************************************************************************/
            /*   Start the loops on wavenumber. This loop will interpolate source terms for scaler Cl	at the discritized nodes			  */
            /*****************************************************************************************************************************/

            if (Variables.itflag != 2)
            {
                if (Variables.twofield == 0)
                {
                    s2 = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sp2 = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sk2 = (double *)malloc((nstep0 + 1) * sizeof(double));
                }
                else if (Variables.twofield == 1)
                {
                    s2adb = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sp2adb = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sk2adb = (double *)malloc((nstep0 + 1) * sizeof(double));

                    s2iso = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sp2iso = (double *)malloc((nstep0 + 1) * sizeof(double));
                    sk2iso = (double *)malloc((nstep0 + 1) * sizeof(double));
                }
            }

            if (Variables.itflag != 0)
            {
                st2 = (double *)malloc((nstep0 + 1) * sizeof(double));
                ste2 = (double *)malloc((nstep0 + 1) * sizeof(double));
                stb2 = (double *)malloc((nstep0 + 1) * sizeof(double));
            }

            if (Variables.itflag != 2)
            {
                if (Variables.twofield == 0)
                {
                    dl2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                }
                if (Variables.twofield == 1)
                {
                    dl2adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl2adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl2adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                    dl2iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl2iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl2iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                }
            }

            //Temporary variables for interpolation
            double ddt, ddt1, ddt2;
            double ddp, ddp1, ddp2;
            double ddk, ddk1, ddk2;

            //Temporary variables for interpolation
            double ddtadb, ddpadb, ddkadb;
            double ddt1adb, ddp1adb, ddk1adb;
            double ddt2adb, ddp2adb, ddk2adb;

            double ddtiso, ddpiso, ddkiso;
            double ddt1iso, ddp1iso, ddk1iso;
            double ddt2iso, ddp2iso, ddk2iso;

            double x, xi, dtau3;
            double h2, a2, b2;

            //Temporary variables for interpolation  Tensor terms
            double ctkj, ctekj, ctbkj, ctckj;

            //Temporary variables for interpolation Tensor terms
            double dtdt, dtde, dtdb;
            double dtdt1, dtde1, dtdb1;
            double dtdt2, dtde2, dtdb2;

            double apowert;

                    
            if (Variables.itflag != 2)
            {

                if (Variables.twofield == 0)
                {
                    dl3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                    dl = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                }
                if (Variables.twofield == 1)
                {
                    dl3adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl3adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl3adb = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                    dl3iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpl3iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkl3iso = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                    dladb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpladb = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkladb = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                    dliso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dpliso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                    dkliso = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                }
                if (Variables.twofield == 0)
                {
                    cl = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cpl = matrix(Variables.l0 + 1, Variables.nn + 1);
                    ccl = matrix(Variables.l0 + 1, Variables.nn + 1);
                    ckkl = matrix(Variables.l0 + 1, Variables.nn + 1);
                    ctkl = matrix(Variables.l0 + 1, Variables.nn + 1);
                }
                if (Variables.twofield == 1)
                {
                    cladb = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cpladb = matrix(Variables.l0 + 1, Variables.nn + 1);
                    ccladb = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cliso = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cpliso = matrix(Variables.l0 + 1, Variables.nn + 1);
                    ccliso = matrix(Variables.l0 + 1, Variables.nn + 1);
                    clcross = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cplcross = matrix(Variables.l0 + 1, Variables.nn + 1);
                    cclcross = matrix(Variables.l0 + 1, Variables.nn + 1);
                }

                for (int j = 1; j < Variables.l0 + 1; j++)
                {
                    for (int in = 1; in < Variables.nn + 1; in++)
                    {
                        if (Variables.twofield == 0)
                        {
                            cl[j][in] = 0.0;
                            cpl[j][in] = 0.0;
                            ccl[j][in] = 0.0;
                            ckkl[j][in] = 0.0;
                            ctkl[j][in] = 0.0;
                        }
                        if (Variables.twofield == 1)
                        {
                            cladb[j][in] = 0.0;
                            cpladb[j][in] = 0.0;
                            ccladb[j][in] = 0.0;
                            cliso[j][in] = 0.0;
                            cpliso[j][in] = 0.0;
                            ccliso[j][in] = 0.0;
                            clcross[j][in] = 0.0;
                            cplcross[j][in] = 0.0;
                            cclcross[j][in] = 0.0;
                        }
                    }
                }
            }

            if (Variables.itflag != 0)
            {
                dtl = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtel = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtbl = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                dtl2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtel2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtbl2 = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                dtl3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtel3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));
                dtbl3 = (double *)malloc((Variables.l0 + 1) * sizeof(double));

                ctl = matrix(Variables.l0 + 1, Variables.nn + 1);
                ctel = matrix(Variables.l0 + 1, Variables.nn + 1);
                ctbl = matrix(Variables.l0 + 1, Variables.nn + 1);
                ctcl = matrix(Variables.l0 + 1, Variables.nn + 1);
            }

            int klo, khi;
            double ho, b0, akt;
            int nstps, nstpt;
            double xlim, xlmax1, xlmax2;
            double xlimmin = 10.0;
            double tmin, tmax, ajl0;
                
            klo = 1;
            khi = 2;
            //             for (int ik = kkmin; ik <= kkmax; ik++)
            //                for (int j = 2; j <= nstep; j++)
            //                {
            //                    printf("%d %d %e %e %e %e %e %e\n", ik, j, d[ik][j], dp[ik][j], dkk[ik][j], dpr[ik][j], dppr[ik][j], dkpr[ik][j]);
            //                    }
            for (int k = kkmin; k <= kkmax; k++)
            {

                akt = ak1[k]; //	Finding position of k in table ak0 to do the interpolation.

                for (klo = 1, khi = 2; ((akt > ak0[klo + 1]) && (klo < nk - 1)); klo++)
                    ; // Check.. There may be some problem.
                khi = klo + 1;

                ho = ak0[khi] - ak0[klo];
                a0 = (ak0[khi] - akt) / ho;
                b0 = (akt - ak0[klo]) / ho;

                nstps = 0;
                nstpt = 0;

                /***********************************************************************************************************************/
                /*     If scalers wanted then this loop will interpolate the source terms at the nodes for tensor Cl						  */
                /***********************************************************************************************************************/
                if (Variables.itflag != 2) // Variables.itflag=2 -> Only Tensor part
                {

                    // Interpolating the source as a function of time for the present wavelength.
                    if (akt < ak10)
                        nstps = nstep - 1;
                    else
                        nstps = n10;

                    if (Variables.twofield == 0)
                    {
                        s2[1] = 0.0;  // Interpolating the source terms
                        sp2[1] = 0.0; // The interpolated source terms will be temporarily stored in these variables
                        sk2[1] = 0.0; //
                    }
                    else if (Variables.twofield == 1)
                    {
                        s2adb[1] = 0.0;  // Interpolating the source terms
                        sp2adb[1] = 0.0; // The interpolated source terms will be temporarily stored in these variables
                        sk2adb[1] = 0.0; //

                        s2iso[1] = 0.0;  // Interpolating the source terms
                        sp2iso[1] = 0.0; // The interpolated source terms will be temporarily stored in these variables
                        sk2iso[1] = 0.0; //
                    }
                    for (int i = 2; i <= nstps; i++) // Spline interpolation of the source terms
                    {
                        if (Variables.twofield == 0)
                        {
                            s2[i] = a0 * d[klo][i] + b0 * d[khi][i] + ((a0 * a0 * a0 - a0) * dpr[klo][i] + (b0 * b0 * b0 - b0) * dpr[khi][i]) * ho * ho / 6.0;
                            sp2[i] = a0 * dp[klo][i] + b0 * dp[khi][i] + ((a0 * a0 * a0 - a0) * dppr[klo][i] + (b0 * b0 * b0 - b0) * dppr[khi][i]) * ho * ho / 6.0;
                            sk2[i] = a0 * dkk[klo][i] + b0 * dkk[khi][i] + ((a0 * a0 * a0 - a0) * dkpr[klo][i] + (b0 * b0 * b0 - b0) * dkpr[khi][i]) * ho * ho / 6.0;
                        }
                        else if (Variables.twofield == 1)
                        {
                            s2adb[i] = a0 * dadb[klo][i] + b0 * dadb[khi][i] + ((a0 * a0 * a0 - a0) * dpradb[klo][i] + (b0 * b0 * b0 - b0) * dpradb[khi][i]) * ho * ho / 6.0;
                            sp2adb[i] = a0 * dpadb[klo][i] + b0 * dpadb[khi][i] + ((a0 * a0 * a0 - a0) * dppradb[klo][i] + (b0 * b0 * b0 - b0) * dppradb[khi][i]) * ho * ho / 6.0;
                            sk2adb[i] = a0 * dkkadb[klo][i] + b0 * dkkadb[khi][i] + ((a0 * a0 * a0 - a0) * dkpradb[klo][i] + (b0 * b0 * b0 - b0) * dkpradb[khi][i]) * ho * ho / 6.0;

                            s2iso[i] = a0 * diso[klo][i] + b0 * diso[khi][i] + ((a0 * a0 * a0 - a0) * dpriso[klo][i] + (b0 * b0 * b0 - b0) * dpriso[khi][i]) * ho * ho / 6.0;
                            sp2iso[i] = a0 * dpiso[klo][i] + b0 * dpiso[khi][i] + ((a0 * a0 * a0 - a0) * dppriso[klo][i] + (b0 * b0 * b0 - b0) * dppriso[khi][i]) * ho * ho / 6.0;
                            sk2iso[i] = a0 * dkkiso[klo][i] + b0 * dkkiso[khi][i] + ((a0 * a0 * a0 - a0) * dkpriso[klo][i] + (b0 * b0 * b0 - b0) * dkpriso[khi][i]) * ho * ho / 6.0;
                        }
                    }

                    if (Variables.twofield == 0)
                    {
                        s2[nstps + 1] = 0.0;
                        sp2[nstps + 1] = 0.0;
                        sk2[nstps + 1] = 0.0;
                    }
                    else if (Variables.twofield == 1)
                    {
                        s2adb[nstps + 1] = 0.0;
                        sp2adb[nstps + 1] = 0.0;
                        sk2adb[nstps + 1] = 0.0;

                        s2iso[nstps + 1] = 0.0;
                        sp2iso[nstps + 1] = 0.0;
                        sk2iso[nstps + 1] = 0.0;
                    }
                }
                /***********************************************************************************************************************/
                /*     If tensors wanted then this loop will interpolate the source terms at the nodes for tensor Cl						  */
                /***********************************************************************************************************************/

                double xf;

                if (Variables.itflag != 0)
                {

                    // Interpolating the tensor source as a function of time for the present wavelength.
                    st2[1] = 0.0;
                    ste2[1] = 0.0;
                    stb2[1] = 0.0; // Initiallize the source terms
                    nstpt = 2;

                    for (int i = 2; i <= nstep; i++)
                    {
                        xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                        if (((akt * GlobalArray.atau0[i]) < stpt) && (xf > 1.0e-8)) // Spline interpolation of the source terms
                        {
                            nstpt = i;
                            st2[i] = a0 * dt[klo][i] + b0 * dt[khi][i] + ((a0 * a0 * a0 - a0) * dtpr[klo][i] + (b0 * b0 * b0 - b0) * dtpr[khi][i]) * ho * ho / 6.0;
                            ste2[i] = a0 * dte[klo][i] + b0 * dte[khi][i] + ((a0 * a0 * a0 - a0) * dtepr[klo][i] + (b0 * b0 * b0 - b0) * dtepr[khi][i]) * ho * ho / 6.0;
                            stb2[i] = a0 * dtb[klo][i] + b0 * dtb[khi][i] + ((a0 * a0 * a0 - a0) * dtbpr[klo][i] + (b0 * b0 * b0 - b0) * dtbpr[khi][i]) * ho * ho / 6.0;
                        }
                        else
                        {
                            st2[i] = 0.0;
                            ste2[i] = 0.0;
                            stb2[i] = 0.0;
                        }
                    }

                    nstpt = (nstpt > n1) ? nstpt : n1;
                    st2[nstpt] = 0.0;
                    ste2[nstpt] = 0.0;
                    stb2[nstpt] = 0.0;
                }

                // (P16) : store k*tau in mxx() for calculating the bessel functions

                double de2;
                /***********************************************************************************************************************/
                /*    Findind the position in the xx table for the x correponding to each timestep												  */
                /***********************************************************************************************************************/
                for (int i = 1; i <= Numericx.maxint(nstps, nstpt) + 2; i++)
                {
                    xf = fabs(akt * (Variables.tau0 - GlobalArray.atau0[i])); // Store k*t values at the grid points

                    if (xf <= 5.0)
                        de2 = 10.0 * xf + 1.0;

                    else if (xf <= 25.0 && xf > 5.0)
                        de2 = (xf - 5.0) * 5.0 + 51.0;

                    else
                        de2 = (xf - 25.0) + 151.0;

                    mxx[i] = (int)(de2);
                    mxx[i] = Numericx.maxint(mxx[i], 1); // Temporarily store k*t values at the grid points
                }

                /*********************************************************************************************************/
                /*     Begin l and  time-loop to integrate scalar perturbations.														*/
                /*     Determining ranges of integration																						*/
                /*********************************************************************************************************/
                // (P17) : Calculating the Cl

                if (Variables.itflag != 2)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        if (Variables.twofield == 0)
                        {
                            dl[j] = 0.0;
                            dl2[j] = 0.0;
                            dl3[j] = 0.0;
                            dpl[j] = 0.0;
                            dpl2[j] = 0.0;
                            dpl3[j] = 0.0;
                            dkl[j] = 0.0;
                            dkl2[j] = 0.0;
                            dkl3[j] = 0.0;
                        }
                        else if (Variables.twofield == 1)
                        {
                            dladb[j] = 0.0;
                            dl2adb[j] = 0.0;
                            dl3adb[j] = 0.0;
                            dpladb[j] = 0.0;
                            dpl2adb[j] = 0.0;
                            dpl3adb[j] = 0.0;
                            dkladb[j] = 0.0;
                            dkl2adb[j] = 0.0;
                            dkl3adb[j] = 0.0;

                            dliso[j] = 0.0;
                            dl2iso[j] = 0.0;
                            dl3iso[j] = 0.0;
                            dpliso[j] = 0.0;
                            dpl2iso[j] = 0.0;
                            dpl3iso[j] = 0.0;
                            dkliso[j] = 0.0;
                            dkl2iso[j] = 0.0;
                            dkl3iso[j] = 0.0;
                        }
                    }
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        xlim = 0.05 * Variables.l[j];
                        xlim = Numericx.max(xlim, xlimmin);
                        xlim = Variables.l[j] - xlim;
                        xlmax1 = 80.0 * Variables.l[j];
                        xlmax2 = Numericx.min(2.0 * Variables.l[j], 1.0 * kmaxjl);

                        tmin = Variables.tau0 - xlmax1 / akt;
                        tmax = Variables.tau0 - xlim / akt;
                        tmax = Numericx.min(Variables.tau0, tmax);
                        tmin = Numericx.max(GlobalArray.atau0[2], tmin);

                        if (tmax < GlobalArray.atau0[2])
                        {
                            goto MOVETOA;
                        }
                        if (Variables.zri == 0.0)
                        {
                            if (tmin < taurend)
                                nstart1 = 2;
                            else
                                nstart1 = n1 + (int)(log(tmin / taurend) / dlntau0);

                            if (tmax < taurend)
                            {
                                nstop1 = n1;
                                nstop1a = n1;
                            }
                            else
                            {
                                nstop1 = n1 + (int)(log(tmax / taurend) / dlntau0);
                                nstop1 = Numericx.minint(nstop1, nstps);

                                if ((akt * GlobalArray.dtau1[nstop1]) > 1.0)
                                    nstop1a = Numericx.maxint(n1, n1 - int(log(akt * dlntau0 * taurend) / dlntau0));
                                else
                                    nstop1a = nstop1;
                            }
                        }
                        else
                        {
                            if (tmin < taurend)
                            {
                                nstart1 = 2;
                                nstart2 = GlobalArray.j2ri1;
                            }
                            else
                            {
                                if (tmin < GlobalArray.atau0[nriend])
                                {
                                    nstart1 = Numericx.minint(n1 + (int)(log((tmin / taurend)) / dlntau0), GlobalArray.j2ri1);
                                    nstart2 = GlobalArray.j2ri1;
                                }
                                else
                                {
                                    nstart1 = nstep + 1;
                                    nstart2 = nriend + (int)(log(tmin / GlobalArray.atau0[nriend] / dlntau0));
                                }
                            }
                            if (tmax < taurend)
                            {
                                nstop1 = n1;
                                nstop2 = 0;
                            }
                            else
                            {
                                if (tmax < GlobalArray.atau0[GlobalArray.j2ri1])
                                {
                                    nstop1 = n1 + (int)(log(tmax / taurend) / dlntau0);
                                    nstop2 = 0;
                                }
                                else
                                {
                                    nstop1 = GlobalArray.j2ri1 - 1;
                                    nstop2 = Numericx.maxint(nriend, nriend + (int)(log(tmax / GlobalArray.atau0[nriend]) / dlntau0));

                                    nstop2 = Numericx.minint(nstop2, nstps);

                                    if ((akt * GlobalArray.dtau1[nstop2]) > 1)
                                        nstop2a = Numericx.maxint(nriend, nriend - (int)(log(akt * dlntau0 * GlobalArray.atau0[nriend]) / dlntau0));
                                    else
                                        nstop2a = nstop2;
                                }
                            
                                nstop1 = Numericx.minint(nstop1, GlobalArray.j2ri1 - 1);

                                if ((akt * GlobalArray.dtau1[nstop1]) > 1)
                                    nstop1a = Numericx.maxint(n1, n1 - (int)(log(akt * dlntau0 * taurend) / dlntau0));
                                else
                                    nstop1a = nstop1;

                            }
                        }

                        /***************************************************************************************************************************/
                        /*    Integration before reionization																													*/
                        /*    Interpolating jls at points where the sources are recorded.																				*/
                        /***************************************************************************************************************************/
                        char c;
                        for (int i = nstart1; i <= nstop1a; i++)
                        {

                            xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                            m2 = mxx[i];
                            h2 = xx[m2 + 1] - xx[m2];
                            a2 = (xx[m2 + 1] - xf) / h2;
                            b2 = (xf - xx[m2]) / h2;
                            ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;

                            if (Variables.twofield == 0)
                            {
                                dl2[j] = dl2[j] + s2[i] * ajl0 * GlobalArray.dtau2[i];
                                dpl2[j] = dpl2[j] + sp2[i] * ajl0 * GlobalArray.dtau2[i];
                                dkl2[j] = dkl2[j] + sk2[i] * ajl0 * GlobalArray.dtau2[i];
                            }

                            if (Variables.twofield == 1)
                            {
                                dl2adb[j] = dl2adb[j] + s2adb[i] * ajl0 * GlobalArray.dtau2[i];
                                dpl2adb[j] = dpl2adb[j] + sp2adb[i] * ajl0 * GlobalArray.dtau2[i];
                                dkl2adb[j] = dkl2adb[j] + sk2adb[i] * ajl0 * GlobalArray.dtau2[i];

                                dl2iso[j] = dl2iso[j] + s2iso[i] * ajl0 * GlobalArray.dtau2[i];
                                dpl2iso[j] = dpl2iso[j] + sp2iso[i] * ajl0 * GlobalArray.dtau2[i];
                                dkl2iso[j] = dkl2iso[j] + sk2iso[i] * ajl0 * GlobalArray.dtau2[i];
                            }
                        }
 
                        /***************************************************************************************************************************/
                        /*     Scalar case.																																			*/
                        /*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
                        /***************************************************************************************************************************/
                        if (Variables.twofield == 0)
                        {
                            for (int i = nstop1a + 1; i <= nstop1; i++)
                            {
                                xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                m2 = mxx[i];
                                dtau3 = GlobalArray.dtau1[i] * fabs(s2[i] / (s2[i + 1] - s2[i] + 1.0e-10));
                                if ((xf < xlmax2) || ((akt * dtau3) < 1.0))
                                {

                                    xi = xf - akt * GlobalArray.dtau1[i];
                                    m1 = mxx[i + 1];
                                    ddt1 = s2[i];
                                    ddt2 = s2[i + 1];
                                    ddp1 = sp2[i];
                                    ddp2 = sp2[i + 1];
                                    ddk1 = sk2[i];
                                    ddk2 = sk2[i + 1];
                                    for (int lx = m1 + 1; lx <= m2; lx++) // convert the source terms from function of (k,t) to the function of (j,l)
                                    {                                     // M.Zaldarriaga's thesis Eq(3.37)
                                        x = xx[lx];
                                        ddt = (ddt1 - ddt2) * (x - xi) / (xf - xi) + ddt2;
                                        ddp = (ddp1 - ddp2) * (x - xi) / (xf - xi) + ddp2;
                                        ddk = (ddk1 - ddk2) * (x - xi) / (xf - xi) + ddk2;
                                        dl3[j] = dl3[j] + ajl[lx][j] * ddt * dxx[lx];
                                        dpl3[j] = dpl3[j] + ajl[lx][j] * ddp * dxx[lx];
                                        dkl3[j] = dkl3[j] + ajl[lx][j] * ddk * dxx[lx];
                                    }
                                }
                            }
                        }
                        else if (Variables.twofield == 1)
                        {
                            for (int i = nstop1a + 1; i <= nstop1; i++)
                            {
                                xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                m2 = mxx[i];
                                dtau3 = GlobalArray.dtau1[i] * fabs(s2adb[i] / (s2adb[i + 1] - s2adb[i] + 1.0e-10));

                                if ((xf < xlmax2) || ((akt * dtau3) < 1.0))
                                {

                                    xi = xf - akt * GlobalArray.dtau1[i];
                                    m1 = mxx[i + 1];

                                    ddt1adb = s2adb[i];
                                    ddt2adb = s2adb[i + 1];
                                    ddp1adb = sp2adb[i];
                                    ddp2adb = sp2adb[i + 1];
                                    ddk1adb = sk2adb[i];
                                    ddk2adb = sk2adb[i + 1];

                                    ddt1iso = s2iso[i];
                                    ddt2iso = s2iso[i + 1];
                                    ddp1iso = sp2iso[i];
                                    ddp2iso = sp2iso[i + 1];
                                    ddk1iso = sk2iso[i];
                                    ddk2iso = sk2iso[i + 1];
                                            	 				
                                    for (int lx = m1 + 1; lx <= m2; lx++) // convert the source terms from function of (k,t) to the function of (j,l)
                                    {                                     // M.Zaldarriaga's thesis Eq(3.37)
                                        x = xx[lx];

                                        ddtadb = (ddt1adb - ddt2adb) * (x - xi) / (xf - xi) + ddt2adb;
                                        ddpadb = (ddp1adb - ddp2adb) * (x - xi) / (xf - xi) + ddp2adb;
                                        ddkadb = (ddk1adb - ddk2adb) * (x - xi) / (xf - xi) + ddk2adb;
                                        dl3adb[j] = dl3adb[j] + ajl[lx][j] * ddt * dxx[lx];
                                        dpl3adb[j] = dpl3adb[j] + ajl[lx][j] * ddpadb * dxx[lx];
                                        dkl3adb[j] = dkl3adb[j] + ajl[lx][j] * ddkadb * dxx[lx];

                                        ddtiso = (ddt1iso - ddt2iso) * (x - xi) / (xf - xi) + ddt2iso;
                                        ddpiso = (ddp1iso - ddp2iso) * (x - xi) / (xf - xi) + ddp2iso;
                                        ddkiso = (ddk1iso - ddk2iso) * (x - xi) / (xf - xi) + ddk2iso;
                                        dl3iso[j] = dl3iso[j] + ajl[lx][j] * ddtiso * dxx[lx];
                                        dpl3iso[j] = dpl3iso[j] + ajl[lx][j] * ddpiso * dxx[lx];
                                        dkl3iso[j] = dkl3iso[j] + ajl[lx][j] * ddkiso * dxx[lx];
                                    }
                                }
                            }
                        }
                        /***************************************************************************************************************************/
                        /*    Integration after reionization																													*/
                        /*		Interpolating jls at points where the sources are recorded.																				*/
                        /***************************************************************************************************************************/
                        if (Variables.zri != 0.0)
                        {
                            if (Variables.twofield == 0)
                            {
                                for (int i = nstart2; i <= nstop2a; i++)
                                {
                                    xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                    m2 = mxx[i];
                                    h2 = xx[m2 + 1] - xx[m2];
                                    a2 = (xx[m2 + 1] - xf) / h2;
                                    b2 = (xf - xx[m2]) / h2;
                                    ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;
                                    dl2[j] = dl2[j] + s2[i] * ajl0 * GlobalArray.dtau2[i];
                                    dpl2[j] = dpl2[j] + sp2[i] * ajl0 * GlobalArray.dtau2[i];
                                    dkl2[j] = dkl2[j] + sk2[i] * ajl0 * GlobalArray.dtau2[i];                            
                                }
                            }

                            else if (Variables.twofield == 1)
                            {
                                for (int i = nstart2; i <= nstop2a; i++)
                                {
                                    xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                    m2 = mxx[i];
                                    h2 = xx[m2 + 1] - xx[m2];
                                    a2 = (xx[m2 + 1] - xf) / h2;
                                    b2 = (xf - xx[m2]) / h2;
                                    ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;
                                    dl2adb[j] = dl2adb[j] + s2adb[i] * ajl0 * GlobalArray.dtau2[i];
                                    dpl2adb[j] = dpl2adb[j] + sp2adb[i] * ajl0 * GlobalArray.dtau2[i];
                                    dkl2adb[j] = dkl2adb[j] + sk2adb[i] * ajl0 * GlobalArray.dtau2[i];
                                }
                                for (int i = nstart2; i <= nstop2a; i++)
                                {
                                    xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                    m2 = mxx[i];
                                    h2 = xx[m2 + 1] - xx[m2];
                                    a2 = (xx[m2 + 1] - xf) / h2;
                                    b2 = (xf - xx[m2]) / h2;
                                    ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;
                                    dl2iso[j] = dl2iso[j] + s2iso[i] * ajl0 * GlobalArray.dtau2[i];
                                    dpl2iso[j] = dpl2iso[j] + sp2iso[i] * ajl0 * GlobalArray.dtau2[i];
                                    dkl2iso[j] = dkl2iso[j] + sk2iso[i] * ajl0 * GlobalArray.dtau2[i];
                                }
                            }

                            /***************************************************************************************************************************/
                            /*     Scaler case.																																			*/
                            /*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
                            /***************************************************************************************************************************/
                            for (int i = nstop2a + 1; i <= nstop2; i++)
                            {
                                xi = xf - akt * GlobalArray.dtau1[i];
                                xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                m2 = mxx[i];
                                if (Variables.twofield == 0) {
                                    dtau3 = GlobalArray.dtau1[i] * fabs(s2[i] / (s2[i + 1] - s2[i] + 1.0e-10));
                                }

                                else if (Variables.twofield == 1) {// Check may be something wrong here..
                                    dtau3 = GlobalArray.dtau1[i] * fabs(s2adb[i] / (s2adb[i + 1] - s2adb[i] + 1.0e-10));
                                }

                                if ((xf < xlmax2) || ((akt * dtau3) < 1.0))
                                {
                                    m1 = mxx[i + 1];

                                    if (Variables.twofield == 0)
                                    {
                                        ddt1 = s2[i];
                                        ddt2 = s2[i + 1];
                                        ddp1 = sp2[i];
                                        ddp2 = sp2[i + 1];
                                        ddk1 = sk2[i];
                                        ddk2 = sk2[i + 1];
                                    }
                                    else if (Variables.twofield == 1)
                                    {
                                        ddt1adb = s2adb[i];
                                        ddt2adb = s2adb[i + 1];
                                        ddp1adb = sp2adb[i];
                                        ddp2adb = sp2adb[i + 1];
                                        ddk1adb = sk2adb[i];
                                        ddk2adb = sk2adb[i + 1];

                                        ddt1iso = s2iso[i];
                                        ddt2iso = s2iso[i + 1];
                                        ddp1iso = sp2iso[i];
                                        ddp2iso = sp2iso[i + 1];
                                        ddk1iso = sk2iso[i];
                                        ddk2iso = sk2iso[i + 1];
                                    }


                                    for (int lx = m1 + 1; lx <= m2; lx++) // convert the source terms from function of (k,t) to the function of (j,l)
                                    {                                     // M.Zaldarriaga's thesis Eq(3.37)
                                        x = xx[lx];

                                        if (Variables.twofield == 0)
                                        {
                                            ddt = (ddt1 - ddt2) * (x - xi) / (xf - xi) + ddt2;
                                            ddp = (ddp1 - ddp2) * (x - xi) / (xf - xi) + ddp2;
                                            ddk = (ddk1 - ddk2) * (x - xi) / (xf - xi) + ddk2;
                                            dl3[j] = dl3[j] + ajl[lx][j] * ddt * dxx[lx];
                                            dpl3[j] = dpl3[j] + ajl[lx][j] * ddp * dxx[lx];
                                            dkl3[j] = dkl3[j] + ajl[lx][j] * ddk * dxx[lx];
                                        }
                                        if (Variables.twofield == 1)
                                        {
                                            ddtadb = (ddt1adb - ddt2adb) * (x - xi) / (xf - xi) + ddt2adb;
                                            ddpadb = (ddp1adb - ddp2adb) * (x - xi) / (xf - xi) + ddp2adb;
                                            ddkadb = (ddk1adb - ddk2adb) * (x - xi) / (xf - xi) + ddk2adb;
                                            dl3adb[j] = dl3adb[j] + ajl[lx][j] * ddtadb * dxx[lx];
                                            dpl3adb[j] = dpl3adb[j] + ajl[lx][j] * ddpadb * dxx[lx];
                                            dkl3adb[j] = dkl3adb[j] + ajl[lx][j] * ddkadb * dxx[lx];

                                            ddtiso = (ddt1iso - ddt2iso) * (x - xi) / (xf - xi) + ddt2iso;
                                            ddpiso = (ddp1iso - ddp2iso) * (x - xi) / (xf - xi) + ddp2iso;
                                            ddkiso = (ddk1iso - ddk2iso) * (x - xi) / (xf - xi) + ddk2iso;
                                            dl3iso[j] = dl3iso[j] + ajl[lx][j] * ddtiso * dxx[lx];
                                            dpl3iso[j] = dpl3iso[j] + ajl[lx][j] * ddpiso * dxx[lx];
                                            dkl3iso[j] = dkl3iso[j] + ajl[lx][j] * ddkiso * dxx[lx];
                                        }
                                    } // LOOP OVER lx ENDS HERE
                                }     // IF END HERE
                            }         // FOR LOOP OVER I ENDS HERE

                        }             // IF OF ZRI ENDS HERE

                        if (Variables.twofield == 0)
                        {
                            dl[j] = dl2[j] + dl3[j] / akt;
                            dpl[j] = dpl2[j] + dpl3[j] / akt;
                            dkl[j] = dkl2[j] + dkl3[j] / akt;
                            //		printf("Hi, %e\n",dl[j]);
                        }
                        if (Variables.twofield == 1)
                        {
                            dladb[j] = dl2adb[j] + dl3adb[j] / akt;
                            dpladb[j] = dpl2adb[j] + dpl3adb[j] / akt;
                            dkladb[j] = dkl2adb[j] + dkl3adb[j] / akt;

                            dliso[j] = dl2iso[j] + dl3iso[j] / akt;
                            dpliso[j] = dpl2iso[j] + dpl3iso[j] / akt;
                            dkliso[j] = dkl2iso[j] + dkl3iso[j] / akt;
                            //printf("%e %e\n",dladb[j],dliso[j]);
                        }

                        //	printf("%d\t %d\t %e\t %e\t %e\n",k,j,dkl2[j],dkl3[j],dkl[j]);
                    }
                    //FOR LOOP OVER J ENDS HERE
                    
                MOVETOA:

                    /***************************************************************************************************************************/
                    /*     Scaler case.																																			*/
                    /*		 Finally calculate Cl by intigrating over K's . ps: l*(l-1) will be done in the next step										*/
                    /***************************************************************************************************************************/

                    //printf("Variables.nn = %d %d %d\n",Variables.nn,Variables.l0,Variables.twofield);
                    double ckj, cpkj, cckj, ckkkj, ctkkj;
                    for (int in = 1; in <= Variables.nn; in++)
                    {
                        if (Variables.twofield == 0)
                        {
                            for (int j = 1; j <= Variables.l0; j++)
                            {
                                Scalar.powersflat(akt, in, &apowers); // apower is the initial power spectrum
                                apowers = apowers / akt;

                                ckj = apowers * dl[j] * dl[j] * dak1[k];
                                cl[j][in] = cl[j][in] + ckj;

                                cpkj = apowers * dpl[j] * dpl[j] * dak1[k]; // Calculating Cl's
                                cpl[j][in] = cpl[j][in] + cpkj;             // Zaldarriaga PhD Thesis 3.25

                                cckj = apowers * dl[j] * dpl[j] * dak1[k];
                                ccl[j][in] = ccl[j][in] + cckj;

                                ckkkj = apowers * dkl[j] * dkl[j] * dak1[k];
                                ckkl[j][in] = ckkl[j][in] + ckkkj;

                                ctkkj = apowers * dl[j] * dkl[j] * dak1[k];
                                ctkl[j][in] = ctkl[j][in] + ctkkj;
                                //printf("%e\t%e\t%e\%e\n",cl[j][in],apowers,dl[j],dak1[k]);
                            }
                        }
                        else if (Variables.twofield == 1)
                        {
                            if (Variables.powerfileflag == 1)
                                Doubleinflation.interpoletepowers(akt);
                            //printf("Hi %d\n",Variables.twofield);
                            for (int j = 1; j <= Variables.l0; j++)
                            {
                                //Others.powersflat(akt,in,&apowers);
                                if (Variables.powerfileflag == 0)
                                    apowers = Doubleinflation.poweradiabatic(akt); // apower is the initial power spectrum
                                else if (Variables.powerfileflag == 1)
                                    apowers = Doubleinflation.powerA;

                                apowers = apowers / akt;

                                ckj = apowers * dladb[j] * dladb[j] * dak1[k];
                                cladb[j][in] = cladb[j][in] + ckj;

                                cpkj = apowers * dpladb[j] * dpladb[j] * dak1[k]; // Calculating Cl's
                                cpladb[j][in] = cpladb[j][in] + cpkj;             // Zaldarriaga PhD Thesis 3.25

                                cckj = apowers * dladb[j] * dpladb[j] * dak1[k];
                                ccladb[j][in] = ccladb[j][in] + cckj;
                            }

                            for (int j = 1; j <= Variables.l0; j++)
                            {
                                if (Variables.powerfileflag == 0)
                                    apowers = Doubleinflation.powerisotherm(akt); // apower is the initial power spectrum
                                else if (Variables.powerfileflag == 1)
                                    apowers = Doubleinflation.powerI;

                                apowers = apowers / akt;

                                ckj = apowers * dliso[j] * dliso[j] * dak1[k];
                                cliso[j][in] = cliso[j][in] + ckj;

                                cpkj = apowers * dpliso[j] * dpliso[j] * dak1[k]; // Calculating Cl's
                                cpliso[j][in] = cpliso[j][in] + cpkj;             // Zaldarriaga PhD Thesis 3.25

                                cckj = apowers * dliso[j] * dpliso[j] * dak1[k];
                                ccliso[j][in] = ccliso[j][in] + cckj;
                            }

                            for (int j = 1; j <= Variables.l0; j++)
                            {
                                if (Variables.powerfileflag == 0)
                                    apowers = Doubleinflation.powercross(akt); // apower is the initial power spectrum
                                else if (Variables.powerfileflag == 1)
                                    apowers = Doubleinflation.powerC;

                                apowers = apowers / akt;

                                ckj = apowers * dladb[j] * dliso[j] * dak1[k];
                                clcross[j][in] = clcross[j][in] + ckj;

                                cpkj = apowers * dpladb[j] * dpliso[j] * dak1[k]; // Calculating Cl's
                                cplcross[j][in] = cplcross[j][in] + cpkj;         // Zaldarriaga PhD Thesis 3.25

                                cckj = apowers * dladb[j] * dpliso[j] * dak1[k];
                                cclcross[j][in] = cclcross[j][in] + cckj;
                                //printf("%e \t%e \t%e\n",cladb[j][in],cliso[j][in],clcross[j][in]);
                            }
                        }
                    }

                } // IF OF Variables.itflag != 2 IE SCALER PART LOOP ENDS HERE
                /***************************************************************************************************************************/
                /*     Begin l and  time-loop to integrate tensor perturbations.																				*/
                /*     Finding the ranges of integration.																												*/
                /***************************************************************************************************************************/
                if (Variables.itflag != 0)
                {

                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        dtl2[j] = 0.0;
                        dtel2[j] = 0.0;
                        dtbl2[j] = 0.0;
                    }

                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        //printf("E123 : %d\n",j);
                        xlim = 0.05 * Variables.l[j];
                        xlim = Numericx.max(xlim, xlimmin);
                        xlim = Variables.l[j] - xlim;
                        xlmax1 = 80.0 * Variables.l[j];
                        //printf("\nE123b : ");
                        tmin = Variables.tau0 - xlmax1 / akt;
                        tmax = Variables.tau0 - xlim / akt;
                        tmax = Numericx.min(Variables.tau0, tmax);
                        //printf("\nE123c : ");
                        tmin = Numericx.max(GlobalArray.atau0[2], tmin);
                        //printf("\nE123d : ");
                        //printf("\n %e %e",tmax,GlobalArray.atau0[2]);
                        if (tmax < GlobalArray.atau0[2])
                            goto MOVETOB;

                        if (Variables.zri == 0.0)
                        {
                            if (tmin < taurend)
                                nstart1 = 2;
                            else
                                nstart1 = n1 + int(log(tmin / taurend) / dlntau0);

                            if (tmax < taurend)
                            {
                                nstop1 = n1;
                                nstop1a = n1;
                            }
                            else
                            {
                                nstop1 = n1 + (int)(log(tmax / taurend) / dlntau0);
                                nstop1 = Numericx.minint(nstop1, nstpt);

                                if ((akt * GlobalArray.dtau1[nstop1]) > 1)
                                {
                                    nstop1a = Numericx.maxint(n1, n1 - (int)(log(akt * dlntau0 * taurend) / dlntau0));
                                }
                                else
                                {
                                    nstop1a = nstop1;
                                }
                            }
                        }
                        else
                        {
                            if (tmin < taurend)
                            {
                                nstart1 = 2;
                                nstart2 = GlobalArray.j2ri1;
                            }
                            else
                            {
                                if (tmin < GlobalArray.atau0[nriend])
                                {
                                    nstart1 = Numericx.minint(n1 + int(log(tmin / taurend) / dlntau0), GlobalArray.j2ri1);
                                    nstart2 = GlobalArray.j2ri1;
                                }
                                else
                                {
                                    nstart1 = nstep + 1;
                                    nstart2 = nriend + int(log(tmin / GlobalArray.atau0[nriend] / dlntau0));
                                }
                            }

                            if (tmax < taurend)
                            {
                                nstop1 = n1;
                                nstop2 = 0;
                            }
                            else
                            {
                                if (tmax < GlobalArray.atau0[GlobalArray.j2ri1])
                                {
                                    nstop1 = n1 + (int)(log(tmax / taurend) / dlntau0);
                                    nstop2 = 0;
                                }
                                else
                                {
                                    nstop1 = GlobalArray.j2ri1 - 1;
                                    nstop2 = Numericx.maxint(nriend, nriend + (int)(log(tmax / GlobalArray.atau0[nriend]) / dlntau0));

                                    nstop2 = Numericx.minint(nstop2, nstpt);

                                    if ((akt * GlobalArray.dtau1[nstop2]) > 1)
                                    {
                                        nstop2a = Numericx.maxint(nriend, nriend - (int)(log(akt * dlntau0 * GlobalArray.atau0[nriend]) / dlntau0));
                                    }
                                    else
                                    {
                                        nstop2a = nstop2;
                                    }
                                }

                                nstop1 = Numericx.minint(nstop1, GlobalArray.j2ri1 - 1);

                                if ((akt * GlobalArray.dtau1[nstop1]) > 1)
                                    nstop1a = Numericx.maxint(n1, n1 - int(log(akt * dlntau0 * taurend) / dlntau0));

                                else
                                    nstop1a = nstop1;
                            }
                        }
                        /***************************************************************************************************************************/
                        /*    Integration before reionization																													*/
                        /*		Interpolating jls at points where the sources are recorded.																				*/
                        /***************************************************************************************************************************/
                        for (int i = nstart1; i <= nstop1a; i++)
                        {
                            xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);

                            if ((j == 1) && (xf < 1.0) && (xf > 1.0e-7))
                                ajl0 = (3.0 * (sin(xf) / xf - cos(xf)) / xf - sin(xf)) / xf;
                            else
                            {
                                m2 = mxx[i];
                                h2 = xx[m2 + 1] - xx[m2];
                                a2 = (xx[m2 + 1] - xf) / h2;
                                b2 = (xf - xx[m2]) / h2;
                                ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;
                            }

                            dtl2[j] = dtl2[j] + st2[i] * ajl0 * GlobalArray.dtau2[i];
                            dtel2[j] = dtel2[j] + ste2[i] * ajl0 * GlobalArray.dtau2[i];
                            dtbl2[j] = dtbl2[j] + stb2[i] * ajl0 * GlobalArray.dtau2[i];
                        }

                        /***************************************************************************************************************************/
                        /*     Tensor case.																																			*/
                        /*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
                        /***************************************************************************************************************************/
                        for (int i = nstop1a + 1; i <= nstop1; i++)
                        {
                            xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                            m2 = mxx[i];
                            xi = xf - akt * GlobalArray.dtau1[i];
                            m1 = mxx[i + 1];
                            dtdt1 = st2[i];
                            dtdt2 = st2[i + 1];
                            dtde1 = ste2[i];
                            dtde2 = ste2[i + 1];
                            dtdb1 = stb2[i];
                            dtdb2 = stb2[i + 1];

                            for (int lx = m1 + 2; lx <= m2 + 1; lx++)
                            {
                                x = xx[lx];
                                dtdt = (dtdt1 - dtdt2) * (x - xi) / (xf - xi) + dtdt2;
                                dtde = (dtde1 - dtde2) * (x - xi) / (xf - xi) + dtde2;
                                dtdb = (dtdb1 - dtdb2) * (x - xi) / (xf - xi) + dtdb2;
                                dtl3[j] = dtl3[j] + ajl[lx][j] * dtdt * dxx[lx];
                                dtel3[j] = dtel3[j] + ajl[lx][j] * dtde * dxx[lx];
                                dtbl3[j] = dtbl3[j] + ajl[lx][j] * dtdb * dxx[lx];
                            }
                        }

                        /***************************************************************************************************************************/
                        /*    Integration after reionization																													*/
                        /*		Interpolating jls at points where the sources are recorded.																				*/
                        /***************************************************************************************************************************/
                        if (Variables.zri != 0.0)
                        {

                            for (int i = nstart2; i <= nstop2a; i++)
                            {
                                xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);

                                if ((j == 1) && (xf < 1.0) && (xf > 1.0e-7))
                                    ajl0 = (3.0 * (sin(xf) / xf - cos(xf)) / xf - sin(xf)) / xf;
                                else
                                {
                                    m2 = mxx[i];
                                    h2 = xx[m2 + 1] - xx[m2];
                                    a2 = (xx[m2 + 1] - xf) / h2;
                                    b2 = (xf - xx[m2]) / h2;
                                    ajl0 = a2 * ajl[m2][j] + b2 * ajl[m2 + 1][j] + ((a2 * a2 * a2 - a2) * ajlpr[m2][j] + (b2 * b2 * b2 - b2) * ajlpr[m2 + 1][j]) * (h2 * h2) / 6.0;
                                }

                                dtl2[j] = dtl2[j] + st2[i] * ajl0 * GlobalArray.dtau2[i];
                                dtel2[j] = dtel2[j] + ste2[i] * ajl0 * GlobalArray.dtau2[i];
                                dtbl2[j] = dtbl2[j] + stb2[i] * ajl0 * GlobalArray.dtau2[i];
                            }

                            /***************************************************************************************************************************/
                            /*     Tensor case.																																			*/
                            /*		 Convert the source terms as a function of l	by multiplying with Jl(x)																*/
                            /***************************************************************************************************************************/
                            for (int i = nstop2a + 1; i <= nstop2; i++)
                            {
                                xi = xf - akt * GlobalArray.dtau1[i];
                                xf = akt * (Variables.tau0 - GlobalArray.atau0[i]);
                                m2 = mxx[i];
                                m1 = mxx[i + 1];
                                dtdt1 = st2[i];
                                dtdt2 = st2[i + 1];
                                dtde1 = ste2[i];
                                dtde2 = ste2[i + 1];
                                dtdb1 = stb2[i];
                                dtdb2 = stb2[i + 1];

                                for (int lx = m1 + 2; lx <= m2 + 1; lx++)
                                {
                                    x = xx[lx];
                                    dtdt = (dtdt1 - dtdt2) * (x - xi) / (xf - xi) + dtdt2;
                                    dtde = (dtde1 - dtde2) * (x - xi) / (xf - xi) + dtde2;
                                    dtdb = (dtdb1 - dtdb2) * (x - xi) / (xf - xi) + dtdb2;
                                    dtl3[j] = dtl3[j] + ajl[lx][j] * dtdt * dxx[lx];
                                    dtel3[j] = dtel3[j] + ajl[lx][j] * dtde * dxx[lx];
                                    dtbl3[j] = dtbl3[j] + ajl[lx][j] * dtdb * dxx[lx];
                                }
                            }
                        }

                        dtl[j] = dtl2[j] + dtl3[j] / akt;
                        dtel[j] = dtel2[j] + dtel3[j] / akt;
                        dtbl[j] = dtbl2[j] + dtbl3[j] / akt;
                    }

                MOVETOB:

                    /***************************************************************************************************************************/
                    /*     Tensor case.																																			*/
                    /*		 Finally calculate Cl by intigrating over K's . ps: l*(l-1) will be done in the next step										*/
                    /***************************************************************************************************************************/

                    for (int in = 1; in <= Variables.nn; in++)
                    {
                        for (int j = 1; j <= Variables.l0; j++)
                        {
                            //printf("\nEXYZ");
                            Tensor.powertflat(akt, in, &apowert); // apower is the initial power spectrum
                            apowert = apowert / akt;
                            //printf("\n%d %d %e %e",in,j,dtl[j],ctl[j][in]);
                            ctkj = apowert * dtl[j] * dtl[j] * dak1[k]; // Calculating Cl's
                            ctl[j][in] = ctl[j][in] + ctkj;             // Zaldarriaga PhD Thesis 3.25

                            ctekj = apowert * dtel[j] * dtel[j] * dak1[k];
                            ctel[j][in] = ctel[j][in] + ctekj;

                            ctbkj = apowert * dtbl[j] * dtbl[j] * dak1[k];
                            ctbl[j][in] = ctbl[j][in] + ctbkj;

                            ctckj = apowert * dtl[j] * dtel[j] * dak1[k];
                            ctcl[j][in] = ctcl[j][in] + ctckj;
                        }
                    }

                } // IF Variables.itflag !=0 ENDS HERE
            }
            if (Variables.itflag != 2)
            {
                // Free the Brightness fluctuation functions
                if (Variables.twofield == 0)
                {
                    free(dl3);
                    free(dpl3);
                    free(dkl3);

                    free(dl);
                    free(dpl);
                    free(dkl);
                }
                if (Variables.twofield == 1)
                {
                    free(dl3adb);
                    free(dpl3adb);
                    free(dkl3adb);

                    free(dl3iso);
                    free(dpl3iso);
                    free(dkl3iso);

                    free(dladb);
                    free(dpladb);
                    free(dkladb);

                    free(dliso);
                    free(dpliso);
                    free(dkliso);
                }
            }

            if (Variables.itflag != 2) // Only scalars are requested
            {
                free(d);
                free(dp);
                free(dkk);
                free(phik);

                free(dpr);
                free(dppr);
                free(dkpr);

                if (Variables.twofield == 1) // Two fields inflation is requested
                {
                    free(dadb);
                    free(dpadb);
                    free(dkkadb);

                    free(diso);
                    free(dpiso);
                    free(dkkiso);

                    free(dpradb);
                    free(dppradb);
                    free(dkpradb);

                    free(dpriso);
                    free(dppriso);
                    free(dkpriso);
                }
            }

            if (Variables.itflag != 0) //If tensors are required
            {
                free(dt);
                free(dte);
                free(dtb);

                free(dtpr);
                free(dtepr);
                free(dtbpr);
            }



        } // END OF IF ICT !=1 .i.e. when Cls are required

        // (P18) : Calculate the Cl's and then calculate the primes at all these points
        double ctnorm, xl[lmax];
        /******************************************************************************************************************************/
        /*  Final calculations for CMB output.																														*/
        /******************************************************************************************************************************/
        for (int j = 1; j <= Variables.l0; j++)
            xl[j] = Variables.l[j];

        /******************************************************************************************************************************/
        /*    Scalar case																																					*/
        /******************************************************************************************************************************/
        if (Variables.itflag != 2)
        {
            if (Variables.twofield == 0)
            {
                for (int in = 1; in <= Variables.nn; in++)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
                        cl[j][in] = 2.0 * cl[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        cpl[j][in] = 2.0 * ctnorm * cpl[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        ccl[j][in] = 2.0 * sqrt(ctnorm) * ccl[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        ckkl[j][in] = 2.0 * ckkl[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        ctkl[j][in] = 2.0 * ctkl[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        //printf("%e\n",cl[j][in]);
                    }
                }
            }

            if (Variables.twofield == 1)
            {
                for (int in = 1; in <= Variables.nn; in++)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
                        cladb[j][in] = 2.0 * cladb[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        cpladb[j][in] = 2.0 * ctnorm * cpladb[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        ccladb[j][in] = 2.0 * sqrt(ctnorm) * ccladb[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                    }
                }

                for (int in = 1; in <= Variables.nn; in++)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
                        cliso[j][in] = 2.0 * cliso[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        cpliso[j][in] = 2.0 * ctnorm * cpliso[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        ccliso[j][in] = 2.0 * sqrt(ctnorm) * ccliso[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                    }
                }

                for (int in = 1; in <= Variables.nn; in++)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
                        clcross[j][in] = 2.0 * clcross[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        cplcross[j][in] = 2.0 * ctnorm * cplcross[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                        cclcross[j][in] = 2.0 * sqrt(ctnorm) * cclcross[j][in] * Variables.l[j] * (Variables.l[j] + 1);
                    }
                }

                for (int in = 1; in <= Variables.nn; in++)
                {
                    for (int j = 1; j <= Variables.l0; j++)
                    {
                        cltotal[j][in] = -(4.0 / 15) * clcross[j][in] + (1.0 / 9) * cladb[j][in] + (4.0 / 25) * cliso[j][in];
                        cpltotal[j][in] = -(4.0 / 15) * cplcross[j][in] + (1.0 / 9) * cpladb[j][in] + (4.0 / 25) * cpliso[j][in];
                        ccltotal[j][in] = -(4.0 / 15) * cclcross[j][in] + (1.0 / 9) * ccladb[j][in] + (4.0 / 25) * ccliso[j][in];

                        // cltotal[j][in] =clcross[j][in]  + 2*cladb[j][in]  + cliso[j][in];
                        // cpltotal[j][in]=cplcross[j][in] + 2*cpladb[j][in] + cpliso[j][in];
                        // ccltotal[j][in]=cclcross[j][in] + 2*ccladb[j][in] + ccliso[j][in];
                    }
                }
            }

            /***************************************************************************************************************************/
            /*   Making the interpolation tables to get other l-values.																						*/
            /***************************************************************************************************************************/

            llo = 1;
            cllo = 1.0e30;
            clhi = 1.0e30;

            if (Variables.twofield == 0)
            {

                // Clpr for interpolating to final l grid
                clpr = matrix(Variables.l0 + 1, Variables.nn + 1);
                cplpr = matrix(Variables.l0 + 1, Variables.nn + 1);
                cclpr = matrix(Variables.l0 + 1, Variables.nn + 1);
                ckklpr = matrix(Variables.l0 + 1, Variables.nn + 1);
                ctklpr = matrix(Variables.l0 + 1, Variables.nn + 1);
                //for (int in = 1; in <= Variables.nn; in++)
                //{
                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	bpcl[jj] = cl[jj][in];
                // 	bpcpl[jj] = cpl[jj][in];
                // 	bpccl[jj] = ccl[jj][in];
                // 	bpckkl[jj] = ckkl[jj][in];
                // 	bpctkl[jj] = ctkl[jj][in];
                // } 
                Numericx.spline2D(xl, cl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, clpr);
                Numericx.spline2D(xl, cpl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cplpr);
                Numericx.spline2D(xl, ccl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cclpr);
                Numericx.spline2D(xl, ckkl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ckklpr);
                Numericx.spline2D(xl, ctkl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ctklpr);

                //for (int jj = 1; jj <= Variables.l0 + 1; jj++)
                //{
                //    printf("%e %e %e\n",xl[jj],cl[jj][1],clpr[jj][1]);
                //}
                // {
                // 	clpr[jj][in] = bpclpr[jj];
                // 	cplpr[jj][in] = bpcplpr[jj];
                // 	cclpr[jj][in] = bpcclpr[jj];
                // 	ckklpr[jj][in] = bpckklpr[jj];
                // 	ctklpr[jj][in] = bpctklpr[jj];
                // }
                //}
                
            }
            
            if (Variables.twofield == 1)
            {
                //Clpr for interpolating to final l grid

                // Two field inflation  Adb perturbation
                clpradb = matrix(Variables.l0 + 1, Variables.nn + 1);
                cplpradb = matrix(Variables.l0 + 1, Variables.nn + 1);
                cclpradb = matrix(Variables.l0 + 1, Variables.nn + 1);

                // Two field inflatiom Isocurvature term
                clpriso = matrix(Variables.l0 + 1, Variables.nn + 1);
                cplpriso = matrix(Variables.l0 + 1, Variables.nn + 1);
                cclpriso = matrix(Variables.l0 + 1, Variables.nn + 1);

                // Two field inflation cross term
                clprcross = matrix(Variables.l0 + 1, Variables.nn + 1);
                cplprcross = matrix(Variables.l0 + 1, Variables.nn + 1);
                cclprcross = matrix(Variables.l0 + 1, Variables.nn + 1);

                // Total Cl = Abd + Iso + 2* Cross
                clprtotal = matrix(Variables.l0 + 1, Variables.nn + 1);
                cplprtotal = matrix(Variables.l0 + 1, Variables.nn + 1);
                cclprtotal = matrix(Variables.l0 + 1, Variables.nn + 1);

                //					for (int in = 1; in <= Variables.nn; in++)
                //					{
                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	bpcl[jj] = cladb[jj][in];
                // 	bpcpl[jj] = cpladb[jj][in];
                // 	bpccl[jj] = ccladb[jj][in];
                // }

                Numericx.spline2D(xl, cladb, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, clpradb);
                Numericx.spline2D(xl, cpladb, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cplpradb);
                Numericx.spline2D(xl, ccladb, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cclpradb);

                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	clpradb[jj][in] = bpclpr[jj];
                // 	cplpradb[jj][in] = bpcplpr[jj];
                // 	cclpradb[jj][in] = bpcclpr[jj];
                // }
                //					}

                //for (int in = 1; in <= Variables.nn; in++)
                //{
                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	bpcl[jj] = cliso[jj][in];
                // 	bpcpl[jj] = cpliso[jj][in];
                // 	bpccl[jj] = ccliso[jj][in];
                // }

                Numericx.spline2D(xl, cliso, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, clpriso);
                Numericx.spline2D(xl, cpliso, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cplpriso);
                Numericx.spline2D(xl, ccliso, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cclpriso);

                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	clpriso[jj][in] = bpclpr[jj];
                // 	cplpriso[jj][in] = bpcplpr[jj];
                // 	cclpriso[jj][in] = bpcclpr[jj];
                // }
                //}

                //for (int in = 1; in <= Variables.nn; in++)
                //{
                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	bpcl[jj] = clcross[jj][in];
                // 	bpcpl[jj] = cplcross[jj][in];
                // 	bpccl[jj] = cclcross[jj][in];
                // }

                Numericx.spline2D(xl, clcross, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, clprcross);
                Numericx.spline2D(xl, cplcross, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cplprcross);
                Numericx.spline2D(xl, cclcross, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cclprcross);

                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	clprcross[jj][in] = bpclpr[jj];
                // 	cplprcross[jj][in] = bpcplpr[jj];
                // 	cclprcross[jj][in] = bpcclpr[jj];
                // }
                //}

                //for (int in = 1; in <= Variables.nn; in++)
                //{
                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	bpcl[jj] = cltotal[jj][in];
                // 	bpcpl[jj] = cpltotal[jj][in];
                // 	bpccl[jj] = ccltotal[jj][in];
                // }

                Numericx.spline2D(xl, cltotal, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, clprtotal);
                Numericx.spline2D(xl, cpltotal, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cplprtotal);
                Numericx.spline2D(xl, ccltotal, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, cclprtotal);

                // for (int jj = 1; jj <= lmax + 1; jj++)
                // {
                // 	clprtotal[jj][in] = bpclpr[jj];
                // 	cplprtotal[jj][in] = bpcplpr[jj];
                // 	cclprtotal[jj][in] = bpcclpr[jj];
                // }
                //}
            }
        }


        /***************************************************************************************************************************/
        /*     Tensor Case																																			*/
        /***************************************************************************************************************************/
        if (Variables.itflag != 0)
        {
            ctl = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctel = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctbl = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctcl = matrix(Variables.l0 + 1, Variables.nn + 1);

            ctlpr = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctelpr = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctblpr = matrix(Variables.l0 + 1, Variables.nn + 1);
            ctclpr = matrix(Variables.l0 + 1, Variables.nn + 1);

            /************************************************************************************************************************/
            /*    Normalization																																		*/
            /************************************************************************************************************************/
            double ctnorm;
            for (int j = 1; j <= Variables.l0; j++)
            {
                ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
                for (int in = 1; in <= Variables.nn; in++)
                {
                    ctl[j][in] = ctnorm * ctl[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
                    ctel[j][in] = ctel[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
                    ctbl[j][in] = ctbl[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
                    ctcl[j][in] = sqrt(ctnorm) * ctcl[j][in] * Variables.l[j] * (Variables.l[j] + 1) / 8.0;
                }
            }

            /************************************************************************************************************************/
            /*		Making the interpolation tables to get other l-values.																				*/
            /************************************************************************************************************************/
            cllo = 1.0e30;
            clhi = 1.0e30;

            // for (int in = 1; in <= Variables.nn; in++)
            // {
            // 	for (int jj = 1; jj <= lmax + 1; jj++)
            // 	{
            // 		bpctl[jj] = ctl[jj][in];
            // 		bpctel[jj] = ctel[jj][in];
            // 		bpctbl[jj] = ctbl[jj][in];
            // 		bpctcl[jj] = ctcl[jj][in];
            // 	}

            Numericx.spline2D(xl, ctl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ctlpr);
            Numericx.spline2D(xl, ctel, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ctelpr);
            Numericx.spline2D(xl, ctbl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ctblpr);
            Numericx.spline2D(xl, ctcl, Variables.l0 + 1, Variables.nn + 1, cllo, clhi, ctclpr);

            // 	for (int jj = 1; jj <= lmax + 1; jj++)
            // 	{
            // 		ctlpr[jj][in] = bpctlpr[jj];
            // 		ctelpr[jj][in] = bpctelpr[jj];
            // 		ctblpr[jj][in] = bpctblpr[jj];
            // 		ctclpr[jj][in] = bpctclpr[jj];
            // 	}
            // }
        } // IF Variables.itflag!=0 ENDS HERE   // Tensor Case

        /************************************************************************************************************************/
        /*     Calculating Cls for every l.																													*/
        /************************************************************************************************************************/

        double clint, cplint, cclint, ckklint, ctklint;
        double clintadb, cplintadb, cclintadb;
        double clintiso, cplintiso, cclintiso;
        double clintcross, cplintcross, cclintcross;
        double clinttotal, cplinttotal, cclinttotal;
        double ctlint, ctelint, ctblint, ctclint;

        int xi, lhi;
        double ho, b0;
        // (P20) : Interpolate the Cls for all the l's
        for (int in = 1; in <= Variables.nn; in++)
        {
        //printf("\nin: %d nn: %d\n", in, Variables.nn);
            llo = 1;
            for (int il = 2; il <= Variables.l[Variables.l0]; il++)
            {
                xi = il;

                if ((xi > xl[llo + 1]) && (llo < Variables.l0))
                    llo = llo + 1;

                lhi = llo + 1;
                ho = xl[lhi] - xl[llo];
                a0 = (xl[lhi] - xi) / ho;
                b0 = (xi - xl[llo]) / ho;
                //printf("\nIn tensor 3192: %d %d %d",in,il, Variables.l[Variables.l0]);
                if (Variables.itflag != 2)
                {
                    if (Variables.twofield == 0)
                    {
                        clint = a0 * cl[llo][in] + b0 * cl[lhi][in] + ((pow(a0, 3) - a0) * clpr[llo][in] + (pow(b0, 3) - b0) * clpr[lhi][in]) * ho * ho / 6.0;
                        cplint = a0 * cpl[llo][in] + b0 * cpl[lhi][in] + ((pow(a0, 3) - a0) * cplpr[llo][in] + (pow(b0, 3) - b0) * cplpr[lhi][in]) * ho * ho / 6.0;
                        cclint = a0 * ccl[llo][in] + b0 * ccl[lhi][in] + ((pow(a0, 3) - a0) * cclpr[llo][in] + (pow(b0, 3) - b0) * cclpr[lhi][in]) * ho * ho / 6.0;
                        ckklint = a0 * ckkl[llo][in] + b0 * ckkl[lhi][in] + ((pow(a0, 3) - a0) * ckklpr[llo][in] + (pow(b0, 3) - b0) * ckklpr[lhi][in]) * ho * ho / 6.0;
                        ctklint = a0 * ctkl[llo][in] + b0 * ctkl[lhi][in] + ((pow(a0, 3) - a0) * ctklpr[llo][in] + (pow(b0, 3) - b0) * ctklpr[lhi][in]) * ho * ho / 6.0;
                        //printf("%d %d %d", llo, lhi, in);
                        //printf("%d %d %e %e %e %e %e\n", in, il, clint,cplint,cclint,ckklint,ctklint);                
                    }
                    if (Variables.twofield == 1)
                    {
                        clintadb = a0 * cladb[llo][in] + b0 * cladb[lhi][in] + ((pow(a0, 3) - a0) * clpradb[llo][in] + (pow(b0, 3) - b0) * clpradb[lhi][in]) * ho * ho / 6.0;
                        cplintadb = a0 * cpladb[llo][in] + b0 * cpladb[lhi][in] + ((pow(a0, 3) - a0) * cplpradb[llo][in] + (pow(b0, 3) - b0) * cplpradb[lhi][in]) * ho * ho / 6.0;
                        cclintadb = a0 * ccladb[llo][in] + b0 * ccladb[lhi][in] + ((pow(a0, 3) - a0) * cclpradb[llo][in] + (pow(b0, 3) - b0) * cclpradb[lhi][in]) * ho * ho / 6.0;

                        clintiso = a0 * cliso[llo][in] + b0 * cliso[lhi][in] + ((pow(a0, 3) - a0) * clpriso[llo][in] + (pow(b0, 3) - b0) * clpriso[lhi][in]) * ho * ho / 6.0;
                        cplintiso = a0 * cpliso[llo][in] + b0 * cpliso[lhi][in] + ((pow(a0, 3) - a0) * cplpriso[llo][in] + (pow(b0, 3) - b0) * cplpriso[lhi][in]) * ho * ho / 6.0;
                        cclintiso = a0 * ccliso[llo][in] + b0 * ccliso[lhi][in] + ((pow(a0, 3) - a0) * cclpriso[llo][in] + (pow(b0, 3) - b0) * cclpriso[lhi][in]) * ho * ho / 6.0;

                        clintcross = a0 * clcross[llo][in] + b0 * clcross[lhi][in] + ((pow(a0, 3) - a0) * clprcross[llo][in] + (pow(b0, 3) - b0) * clprcross[lhi][in]) * ho * ho / 6.0;
                        cplintcross = a0 * cplcross[llo][in] + b0 * cplcross[lhi][in] + ((pow(a0, 3) - a0) * cplprcross[llo][in] + (pow(b0, 3) - b0) * cplprcross[lhi][in]) * ho * ho / 6.0;
                        cclintcross = a0 * cclcross[llo][in] + b0 * cclcross[lhi][in] + ((pow(a0, 3) - a0) * cclprcross[llo][in] + (pow(b0, 3) - b0) * cclprcross[lhi][in]) * ho * ho / 6.0;

                        clinttotal = a0 * cltotal[llo][in] + b0 * cltotal[lhi][in] + ((pow(a0, 3) - a0) * clprtotal[llo][in] + (pow(b0, 3) - b0) * clprtotal[lhi][in]) * ho * ho / 6.0;
                        cplinttotal = a0 * cpltotal[llo][in] + b0 * cpltotal[lhi][in] + ((pow(a0, 3) - a0) * cplprtotal[llo][in] + (pow(b0, 3) - b0) * cplprtotal[lhi][in]) * ho * ho / 6.0;
                        cclinttotal = a0 * ccltotal[llo][in] + b0 * ccltotal[lhi][in] + ((pow(a0, 3) - a0) * cclprtotal[llo][in] + (pow(b0, 3) - b0) * cclprtotal[lhi][in]) * ho * ho / 6.0;
                    }
                }

                if (Variables.itflag != 0)
                {
                    ctlint = a0 * ctl[llo][in] + b0 * ctl[lhi][in] + ((pow(a0, 3) - a0) * ctlpr[llo][in] + (pow(b0, 3) - b0) * ctlpr[lhi][in]) * ho * ho / 6.0;
                    ctelint = a0 * ctel[llo][in] + b0 * ctel[lhi][in] + ((pow(a0, 3) - a0) * ctelpr[llo][in] + (pow(b0, 3) - b0) * ctelpr[lhi][in]) * ho * ho / 6.0;
                    ctblint = a0 * ctbl[llo][in] + b0 * ctbl[lhi][in] + ((pow(a0, 3) - a0) * ctblpr[llo][in] + (pow(b0, 3) - b0) * ctblpr[lhi][in]) * ho * ho / 6.0;
                    ctclint = a0 * ctcl[llo][in] + b0 * ctcl[lhi][in] + ((pow(a0, 3) - a0) * ctclpr[llo][in] + (pow(b0, 3) - b0) * ctclpr[lhi][in]) * ho * ho / 6.0;

                }




                if (Variables.itflag != 2)
                {
                    if (Variables.twofield == 0)
                    {
                        GlobalArray.clts[il][in] = clint / cl[1][in];
                        GlobalArray.cles[il][in] = cplint / cl[1][in];
                        GlobalArray.clcs[il][in] = cclint / cl[1][in];
                        GlobalArray.clkk[il][in] = ckklint / cl[1][in];
                        GlobalArray.cltk[il][in] = ctklint / cl[1][in];
                    }
                    if (Variables.twofield == 1)
                    {
                        /*
		         	cltsadb[il][in]=clintadb/cladb[1][in];
   	           	clesadb[il][in]=cplintadb/cladb[1][in];
      	        	clcsadb[il][in]=cclintadb/cladb[1][in];

		         	cltsiso[il][in]=clintiso/cliso[1][in];
   	           	clesiso[il][in]=cplintiso/cliso[1][in];
      	        	clcsiso[il][in]=cclintiso/cliso[1][in];

		         	cltscross[il][in]=clintcross/clcross[1][in];
   	           	clescross[il][in]=cplintcross/clcross[1][in];
      	        	clcscross[il][in]=cclintcross/clcross[1][in];
						*/
                        GlobalArray.cltsadb[il][in] = clintadb / cltotal[1][in];
                        GlobalArray.clesadb[il][in] = cplintadb / cltotal[1][in];
                        GlobalArray.clcsadb[il][in] = cclintadb / cltotal[1][in];

                        GlobalArray.cltsiso[il][in] = clintiso / cltotal[1][in];
                        GlobalArray.clesiso[il][in] = cplintiso / cltotal[1][in];
                        GlobalArray.clcsiso[il][in] = cclintiso / cltotal[1][in];

                        GlobalArray.cltscross[il][in] = clintcross / cltotal[1][in];
                        GlobalArray.clescross[il][in] = cplintcross / cltotal[1][in];
                        GlobalArray.clcscross[il][in] = cclintcross / cltotal[1][in];

                        GlobalArray.clts[il][in] = clinttotal / cltotal[1][in];
                        GlobalArray.cles[il][in] = cplinttotal / cltotal[1][in];
                        GlobalArray.clcs[il][in] = cclinttotal / cltotal[1][in];
                    }
                }

                if (Variables.itflag != 0)
                {
                    GlobalArray.cltt[il][in] = ctlint / ctl[1][in];
                    GlobalArray.clet[il][in] = ctelint / ctl[1][in];
                    GlobalArray.clbt[il][in] = ctblint / ctl[1][in];
                    GlobalArray.clct[il][in] = ctclint / ctl[1][in];
                }
            }

            if (Variables.itflag != 2)
            {
                if (Variables.twofield == 0)
                    GlobalArray.clts[1][in] = cl[1][in];
                else if (Variables.twofield == 1)
                {
                    GlobalArray.cltsadb[1][in] = cltotal[1][in];
                    GlobalArray.cltsiso[1][in] = cltotal[1][in];
                    GlobalArray.cltscross[1][in] = cltotal[1][in];
                    GlobalArray.clts[1][in] = cltotal[1][in];
                }
            }

            if (Variables.itflag != 0)
                GlobalArray.cltt[1][in] = ctl[1][in];
        }

        
        printf("\nERROR CMBFLAT #3257");
        
        if (Variables.itflag != 2)
        {
            if (Variables.twofield == 0)
            {
                free(clpr);
                free(cplpr);
                free(cclpr);
                free(ckklpr);
                free(ctklpr);
            }
            if (Variables.twofield == 0)
            {
                free(clpradb);
                free(cplpradb);
                free(cclpradb);

                // Two field inflatiom Isocurvature term
                free(clpriso);
                free(cplpriso);
                free(cclpriso);

                // Two field inflation cross term
                free(clprcross);
                free(cplprcross);
                free(cclprcross);

                // Total Cl = Abd + Iso + 2* Cross
                free(clprtotal);
                free(cplprtotal);
                free(cclprtotal);
            }
            
                if (Variables.twofield == 0)
                {
                    free(cl);
                    free(cpl);
                    free(ccl);
                    free(ckkl);
                    free(ctkl);
                }
                if (Variables.twofield == 1)
                {
                    free(cladb);
                    free(cpladb);
                    free(ccladb);
                    free(cliso);
                    free(cpliso);
                    free(ccliso);
                    free(clcross);
                    free(cplcross);
                    free(cclcross);
                }           
        }
        printf("\nERROR CMBFLAT #3313");
        
        if (Variables.itflag != 0)
        {
            free(ctl);
            free(ctel);
            free(ctbl);
            free(ctcl);

            free(ctlpr);
            free(ctelpr);
            free(ctblpr);
            free(ctclpr);
            }
        
        if (printperturbations == 1)
        {
            fclose(fp[1]);
            fclose(fp[2]);
            fclose(fp[3]);
            fclose(fp[4]);
            fclose(fp[5]);
            fclose(fp[6]);
            fclose(fp[7]);
            fclose(fp[8]);
            fclose(fp[9]);
        }

        printf("taumin: %e\n", taumin);
        printf("taurmax: %e\n", Variables.taurmax);
        printf("taurend: %e\n", taurend);
        printf("taurist: %e\n", Variables.taurist);
        printf("tauristp: %e\n", Variables.tauristp);
        printf("taumax: %e\n", taumax);
        printf("taupresent: %e\n", Variables.tau0);

        FILE *fpclts;
        fpclts = fopen("clts.d", "w");
        for (int i = 0; i < Variables.lmo + 1; i++)
        {
            fprintf(fpclts, "%e \n", GlobalArray.clts[i][1]);
        }
        fclose(fpclts);
        
        printf("\nERROR CMBFLAT #3356");
    }

	/*********************************************************************************/
	/* Normalization is requiredfor matching the data with the observed data set     */
	/* Otherwise different theory will give different amplitude of the spectrum      */
	/* and it would be difficult to find out the accuracy of the theory              */
	/*                                                                               */
	/* This part is based on the paper: 													      */
	/*      The Four-Year COBE Normalization and Large-Scale Structure			      */
	/*      Emory F. Bunn  and  Martin White												      */
	/*********************************************************************************/

	void COBEnormalize()
	{

		const double fourpi = 4.0 * 3.14159265;
		double xlog10, h, xlnh, hc, curv, r;
		double c10, d1, d2, d3, d4, d5, d6, d7, x1, x2, x3, x4, x5, x6, x7;
		double sy, s, sx, sxy, sxx, delt, d1pr, d1ppr;

		double d2norm; // Needed for lensing calculation

		xlog10 = log(10.0);
		h = Variables.h0 / 100.0;
		xlnh = log(h);

		// Curvature radius
		if (fabs(Variables.omegak) > 1.0e-3)
		{
			hc = 2.998e5 / Variables.h0;			// h*c. c = 2.338*e5 km/s. h = km/sec/pc
			curv = -Variables.omegak / (hc * hc); // Gaussian curvature
			r = 1.0 / sqrt(fabs(curv)); // Radious of curvature
		}

		// Ensuring scalar to tensor ratio
		if (Variables.itflag == 1)
			for (int in = 1; in <= Variables.nn; in++)
				for (int il = 2; il <= Variables.l[Variables.l0]; il++)
				{
					printf("rat =  %e\n", Variables.rat[il]);
					GlobalArray.cltt[il][in] = GlobalArray.cltt[il][in] * Variables.rat[in];
					GlobalArray.clet[il][in] = GlobalArray.clet[il][in] * Variables.rat[in];
					GlobalArray.clbt[il][in] = GlobalArray.clbt[il][in] * Variables.rat[in];
					GlobalArray.clct[il][in] = GlobalArray.clct[il][in] * Variables.rat[in];
					GlobalArray.clkk[il][in] = GlobalArray.clkk[il][in] * Variables.rat[in];
					GlobalArray.cltk[il][in] = GlobalArray.cltk[il][in] * Variables.rat[in];
				}

		// COBE normalization
		// fit the spectrum to a quadratic around C_10 with equal weights in log_l

		for (int in = 1; in <= Variables.nn; in++)
		{
			c10 = GlobalArray.clts[10][in] + GlobalArray.cltt[Variables.l[2]][in];
			printf("cl10 = %e cltt[10][in] = %e\n", c10, GlobalArray.clts[10][in]);
			d1 = (GlobalArray.clts[Variables.l[2]][in] + GlobalArray.cltt[Variables.l[2]][in]) / c10 - 1.0;
			d2 = (GlobalArray.clts[Variables.l[3]][in] + GlobalArray.cltt[Variables.l[3]][in]) / c10 - 1.0;
			d3 = (GlobalArray.clts[Variables.l[5]][in] + GlobalArray.cltt[Variables.l[5]][in]) / c10 - 1.0;
			d4 = (GlobalArray.clts[Variables.l[7]][in] + GlobalArray.cltt[Variables.l[7]][in]) / c10 - 1.0;
			d5 = (GlobalArray.clts[Variables.l[10]][in] + GlobalArray.cltt[Variables.l[10]][in]) / c10 - 1.0;
			d6 = (GlobalArray.clts[Variables.l[11]][in] + GlobalArray.cltt[Variables.l[11]][in]) / c10 - 1.0;
			d7 = (GlobalArray.clts[Variables.l[12]][in] + GlobalArray.cltt[Variables.l[12]][in]) / c10 - 1.0;

			x1 = log(1.0 * Variables.l[2]) / xlog10 - 1.0;
			x2 = log(1.0 * Variables.l[3]) / xlog10 - 1.0;
			x3 = log(1.0 * Variables.l[5]) / xlog10 - 1.0;
			x4 = log(1.0 * Variables.l[7]) / xlog10 - 1.0;
			x5 = log(1.0 * Variables.l[10]) / xlog10 - 1.0;
			x6 = log(1.0 * Variables.l[11]) / xlog10 - 1.0;
			x7 = log(1.0 * Variables.l[12]) / xlog10 - 1.0;

			// Quadratic Least Square fit
			sy = x1 * d1 + x2 * d2 + x3 * d3 + x4 * d4 + x5 * d5 + x6 * d6 + x7 * d7;
			s = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 + x6 * x6 + x7 * x7;
			sx = x1 * x1 * x1 + x2 * x2 * x2 + x3 * x3 * x3 + x4 * x4 * x4 + x5 * x5 * x5 + x6 * x6 * x6 + x7 * x7 * x7;
			sxy = x1 * x1 * d1 + x2 * x2 * d2 + x3 * x3 * d3 + x4 * x4 * d4 + x5 * x5 * d5 + x6 * x6 * d6 + x7 * x7 * d7;
			sxx = pow(x1, 4) + pow(x2, 4) + pow(x3, 4) + pow(x4, 4) + pow(x5, 4) + pow(x6, 4) + pow(x7, 4);

			delt = s * sxx - sx * sx;
			d1pr = (sxx * sy - sx * sxy) / delt;
			d1ppr = 2.0 * (s * sxy - sx * sy) / delt;

			printf("%e %e %e\n", delt, d1pr, d1ppr);

			// Bunn and White fitting formula
			c10 = (0.64575 + 0.02282 * d1pr + 0.01391 * d1pr * d1pr - 0.01819 * d1ppr - 0.00646 * d1pr * d1ppr + 0.00103 * d1ppr * d1ppr) / c10;

			// Calculations for lensing

			d2norm = c10 * (1.1e-9) / GlobalArray.clts[1][in] / (2 * 3.14159265);

			if (Variables.ict == 2)
			{
				for (int itf = 1; itf <= Variables.ntf; itf++)
				{
					d2norm = d2norm * exp(xlnh * 4.0);
					Variables.anorm[itf][in] = d2norm / ((4 * 3.14159265));
				}
			}

			// C_l normalization; output l(l+1)C_l/twopi

			c10 = (c10 * 2.2e-9) / fourpi;
			printf("cl10 = %e \n", c10);
			for (int il = 2; il <= Variables.l[Variables.l0]; il++)
			{
				printf("clts = %e \n", GlobalArray.clts[il][in]);
				GlobalArray.clts[il][in] = GlobalArray.clts[il][in] * c10;
				GlobalArray.cles[il][in] = GlobalArray.cles[il][in] * c10;
				GlobalArray.clcs[il][in] = GlobalArray.clcs[il][in] * c10;
				GlobalArray.cltt[il][in] = GlobalArray.cltt[il][in] * c10;
				GlobalArray.clet[il][in] = GlobalArray.clet[il][in] * c10;
				GlobalArray.clbt[il][in] = GlobalArray.clbt[il][in] * c10;
				GlobalArray.clct[il][in] = GlobalArray.clct[il][in] * c10;
				GlobalArray.clkk[il][in] = GlobalArray.clkk[il][in] * c10;
				GlobalArray.cltk[il][in] = GlobalArray.cltk[il][in] * c10;
			}
		}

     
 	}


};

#endif /* FUNCTIONS_H_ */