#include<stdio.h>
#include<math.h>

#include "recombination.h"
#include "variables.h"
#include "cmb.h"


int main() {
    CMB cmb;
    recombination Recombination;
    
    cmb.setparam();
    // print the reionization fraction at redshift z = 1370 (a = 1/(z+1))
    printf("reionization fraction x_H0 = %f", Recombination.ionsaha(3740, .00072939));
}