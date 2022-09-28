#include<stdio.h>
#include<math.h>

#include "recombination.h"
#include "variables.h"
#include "others.h"
#include "reion.h"
#include "numericx.h"
#include "cmb.h"


int main() {
    CMB cmb;
    others Others;
    neutrino Neutrino;
    recombination Recombination;
    reion Reion;
    numericx Numericx;

    cmb.OmegaB = 0.05;
    cmb.OmegaC = 0.25;
    cmb.OmegaDE = 0.7;
    cmb.OmegaNmassive = 0.0;
    cmb.H0 = 67.9;
    cmb.Tcmb = 2.7254;
    cmb.nNeutrinoMassive = 0.0;
    
    cmb.ReionizationFlag = 1;   // Re-ionization Optical depth specified
    cmb.OpticalDepthToReionization = 0.08;
    cmb.RecombinationType = 2;  // Saha recombination

    cmb.setparam();

    /* calculate the neutrino density if massive neutrinos are present */
    if(cmb.OmegaNmassive > 0.0 && cmb.nNeutrinoMassive > 1)
    {
    	Neutrino.initnul();  
    }
    if (Variables.optdlss > 0.0)
        Reion.reiopar(1.0); // It will calculate zri and zristp
    
    if (Variables.zri != 0.0)
    {
    Variables.taurist = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / (1 + Variables.zri), Initval.tol);	 // Conformal time where the reionozation starts
	Variables.tauristp = Numericx.rombint(Dtauda, 1.0e-8, 1.0 / (1 + Variables.zristp), Initval.tol); // Conformal tile where the reionization stops
    }
    else
    {
	Variables.taurist = Variables.tau0;
	Variables.tauristp = Variables.tau0;
	}
        
    double tau_present = Others.calcconftime(1.0e-8);

    int step_size[1];
    double tau_rend[1];
    Variables.tau0 = tau_present;
    printf("I think its working");
    int result = Recombination.finithermo(0.1, tau_present, 0.0025, step_size, tau_rend);
}