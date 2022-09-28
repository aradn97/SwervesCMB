#include<stdio.h>
#include<math.h>

#include "variables.h"
#include "neutrino.h"
#include "cmb.h"

int main() {
    neutrino Neutrino;
    CMB cmb;

    double rhonu[2];
    double rrnu[2];
    
    cmb.OmegaB = 0.05;
    cmb.OmegaC = 0.25;
    cmb.OmegaDE = 0.7-0.0054;
    // mass of non relativistic neutrinos
    cmb.OmegaNmassive = 0.0054;
    cmb.H0 = 67.9;
    cmb.Tcmb = 2.7254;
    cmb.nNeutrinoMassive = 3.0;
    

    cmb.setparam();

    Neutrino.initnul();

    Neutrino.nu1(0.005, rrnu);
    Neutrino.ninul(0.005, rhonu);

    printf("%lf %lf\n", rrnu[0], rrnu[1]);
    printf("%lf %lf\n", rhonu[0], rhonu[1]);
}