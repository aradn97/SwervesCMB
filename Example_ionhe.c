#include<stdio.h>
#include<math.h>

#include "variables.h"
#include "neutrino.h"
#include "recombination.h"
#include "cmb.h"


int main() {
    CMB cmb;
    recombination Recombination;
    
    cmb.setparam();
    double ionHe[2]={0.0,0.0};
    
    Recombination.ionhe(15796.0,0.00017238,1.0,ionHe);
    
    printf("%lf %lf\n", ionHe[0], ionHe[1]);
}