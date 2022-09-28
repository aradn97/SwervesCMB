#include<stdio.h>
#include<math.h>

#include "variables.h"
#include "neutrino.h"
#include "cmb.h"

int main() {
    neutrino Neutrino;
    CMB cmb;

    double rhonu[2];

    FILE *fp;
    /* creating a file to write the results to */
    fp =fopen("NeutrinoMass_Omegan_0054.d","w");
    
    cmb.OmegaB = 0.05;
    cmb.OmegaC = 0.25;
    cmb.OmegaDE = 0.7-0.0054;
    // mass of non relativistic neutrinos
    cmb.OmegaNmassive = 0.0054;
    cmb.H0 = 67.9;
    cmb.Tcmb = 2.7254;
    cmb.nNeutrinoMassive = 3.0;
    

    cmb.setparam();
    for (int i = 10000; i >= 0; i--)
    {   
        // creating a range of scale factor from 1e-5 to 1
        double a = pow(10.0, -5.0*i/10000);
        Neutrino.ninul(a, rhonu);
        fprintf(fp, "%e %e %e\n", a, rhonu[0], rhonu[1]);
    }
    fclose(fp);
}