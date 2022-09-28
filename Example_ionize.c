#include<stdio.h>
#include<math.h>

#include "recombination.h"
#include "variables.h"
#include "others.h"
#include "cmb.h"


int main() {
    CMB cmb;
    others Others;
    recombination Recombination;
    
    cmb.setparam();
    // print the reionization fraction at redshift z = 1370 (a = 1/(z+1))
    double adot = 1/Others.dtauda(.00072939);
    printf("reionization fraction x_e = %f", Recombination.ionize(3740, .00072939, adot, 1, 0.5));
}