#ifndef _CMB_H_
#define _CMB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "initval.h"
#include "variables.h"
#include "neutrino.h"

class CMB
{
    public:

    double OmegaB = 0.0486;
    double OmegaC = 0.2589;
    double OmegaDE = 0.6925;
    double OmegaNmassive = 0.0;
    double H0 = 67.74;
    double Tcmb = 2.7254;
    double YHe = 0.238;
    double nNeutrinoMassLess = 3.04;
    double nNeutrinoMassive = 0;
    double gsMassive = 10.75;
    double w_DE = -1;

    double k_mode = 0.001;                           // k mode

    double ns =0.97;
    double alphans = 0.0;
    double nt = ns - 1.0;
    double alphant = 0.0;
    double Tensor2ScalarRatio = 7.0*(1-ns);


    /*******************************************************/    
    // Reionization variables
    /*******************************************************/        
    double ReionizationFraction = 1.0;              // If Optical depth is an input
    double ReionizationRedShift = 11.3;             // If Reionization Redshift is an input    
    double OpticalDepthToReionization = 0.08;       // Optical Depth to Reionization
            
    /*******************************************************/
    // FRW / DGP(5D model)
    // Default model : FRW    
    // Set DGP = 1 for 5D model
    /*******************************************************/    
    unsigned int DGP = 0;
    
    /*******************************************************/
    // 1 : perturbation
    // 3 : without perturbation
    // 2 : table with a,w with perturbation
    // 4 : table with a,w without perturbation
    // Default Value : 3
    /*******************************************************/    
    unsigned int DE = 3;  			    

    /*******************************************************/
    // 0 : no reionization
    // 1 : specified optical depth to lss(xe=1)
    // 2 : specified redshift and xe
    // Default Value : 1
    /*******************************************************/        
    unsigned int ReionizationFlag = 1;

    /*******************************************************/
   	// 0 : Peebles recombination 
	// 1 : Recfast
	// 2 : Saha
    /*******************************************************/
	unsigned int RecombinationType = 0;

    /*******************************************************/
    // 1 : Adiabatic Initial Condition
    // 2 : CDM Isocurvature initial condition
    // 3 : Baryon Isocurvature initial condition
    // 4 : Isocurvature Seed initial condition
    // 5 : Neutrino Density initial condition
    // 6 : Neutrino Velocity initial condition
    // Default Value : 1
    /*******************************************************/        
    unsigned int InitialCondition = 1;


    /*******************************************************/
    // 0 : For Scalar Alone
    // 1 : For Scalar + Tensor
    // 2 : For Tensor Aloone
    // 3 : For scalar (k<k*)
    // 4 : For Scalar (k>k*)
    // Default Value : 1
    /*******************************************************/        
	unsigned int PerturbationType = 1;

    /*******************************************************/
    // 0: CMB
    // 1: Transfer Functions
    // 2: CMB + Transfer Functions
    /*******************************************************/
    unsigned int CMB_TF_Both = 2;

    /*******************************************************/
    // 0: Single Field Inflation
    // 1: Two Field Inflation
    /*******************************************************/
    unsigned int CMB_Inflation = 0;        
        
    /*******************************************************/
	// Number of multipoles in photon and neutrino perturbation 
    /*******************************************************/	
    unsigned int lmaxg = lmax0;					//Number of photon equations for each polarizations
    unsigned int lmaxnr = lmaxnr0;				 //Number of relativistic neutrino equations 
    unsigned int nqmax = nqmax0;				 //Number of q's at which massive neutrino equations will be written
    unsigned int lmaxnu = lmaxnu0;				 //Number of massive neutrino equations for each q. 	

    unsigned int NumberOfPk = 1;

   // unsigned int Lmax = 1500;


    void setparam()
    {
        Variables.omegab = OmegaB;
        Variables.omegac = OmegaC;
        Variables.omegav = OmegaDE;
        Variables.omegan = OmegaNmassive;


		Variables.lmaxg = lmaxg;				//Number of photon equations for each polarizations
		Variables.lmaxnr = lmaxnr;				 //Number of relativistic neutrino equations 
		Variables.nqmax = nqmax;				 //Number of q's at which massive neutrino equations will be written
		Variables.lmaxnu = lmaxnu;				 //Number of massive neutrino equations for each q. 	




        Variables.dimflag = 0;

        Variables.VPCflag = 0;
        
        if(DGP !=0)
        {
            Variables.dimflag = 1;
        }
        
        Variables.h0 = H0;
        Variables.tcmb = Tcmb;
        Variables.yhe = YHe;
        
        Variables.annur = nNeutrinoMassLess;
        Variables.annunr = nNeutrinoMassive;
        Variables.gsnunr = gsMassive;
        Variables.wdyn = w_DE;
        Variables.ndyn = DE;
        
        Variables.ict = CMB_TF_Both;
        Variables.ntf = 1;      // Number of Transfer Functions required 
        Variables.ztf[1] = 0.0;  // redshift at which the transfer function is requested
        
        if(Variables.ict != 0)  //Both CMB and Transfer functions are requested
        {
            Variables.akmaxt = 5.0;
            Variables.nlnkt = 5;
        }
        
        
        Variables.twofield = CMB_Inflation;  //Inflation Model : Single Field / Two Field     
        Variables.nn = NumberOfPk;                    // Number of P(k) for which it will calculate           

        Variables.riflag = ReionizationFlag;
        if(ReionizationFlag == 1)
        {
            Variables.optdlss = OpticalDepthToReionization;
            Variables.zri = 0.0;
            Variables.rif = 1.0;
        }
        else if(ReionizationFlag == 2)
        {
            Variables.optdlss = 0.0;
            Variables.zri = ReionizationRedShift;
            Variables.rif = ReionizationFraction;
        }
        else
        {
            Variables.optdlss = 0.0;
            Variables.zri = 0.0;
            Variables.rif = 0.0;
        }

		Variables.rcflag = RecombinationType;

		Variables.ak = k_mode;
		Variables.ak2 = Variables.ak*Variables.ak;
		Variables.initfl = InitialCondition;

		Variables.itflag = PerturbationType;
                

        Variables.lensflag = 1;
                

        Variables.setparams(); 


		if ((Variables.annunr != 0) || (Variables.omegan != 0.0))
		{
			neutrino Neutrino;
			Neutrino.initnul();
		}

			for(int i=1;i<=Variables.nn;i++)
			{
				Variables.an[i] = ns;
                Variables.alphans[i] = alphans;

				Variables.ant[i] = nt;
                Variables.alphant[i] = alphant;

                Variables.rat[i] = Tensor2ScalarRatio; 
			}

            GlobalArray.allocate(Variables.lmo,Variables.nn,Variables.twofield);
			}	
    };

    #endif