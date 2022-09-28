#include<stdio.h>
#include<math.h>

#include "variables.h"
#include "others.h"
#include "neutrino.h"
#include "cmb.h"


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
    
    FILE *fptr;
    fptr = fopen("soundspeed.d","w");
     double interval = 0.008;

    for (int i = 0; i <= 1000; i++)
    {   
        double a = pow(10, -(interval*i));
        /* speed of sound in the third column */
        fprintf(fptr, "%e %e %e\n", a, Others.dsoundda(a), Others.dsoundda(a)/Others.dtauda(a));
    }
    fclose(fptr);
}