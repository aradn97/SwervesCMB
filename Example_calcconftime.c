#include<stdio.h>
#include<math.h>

#include "variables.h"
#include "others.h"
#include "neutrino.h"


int main() {
    CMB cmb;
    others Others;
    neutrino Neutrino;

    cmb.OmegaB = 0.05;
    cmb.OmegaC = 0.25;
    cmb.OmegaDE = 0.7;
    cmb.OmegaNmassive = 0.0;
    cmb.H0 = 67.9;
    cmb.Tcmb = 2.7254;
    cmb.nNeutrinoMassive = 0.0;

    cmb.setparam(); 

    /* calculate the neutrino density if massive neutrinos are present */
    if(cmb.OmegaNmassive > 0.0 && cmb.nNeutrinoMassive > 1)
    {
    	Neutrino.initnul();  
    }

    printf("%e",Others.calcconftime(0.5));

}