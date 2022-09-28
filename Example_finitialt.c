#include<stdio.h>
#include<math.h>

#include "others.h"
#include "tensor.h"
#include "cmb.h"

int main() {
    CMB cmb;
    tensor Tensor;
    others Others;

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
    cmb.InitialCondition = 1; // Adiabatic perturbation

	cmb.lmaxg = 12;				//Number of Photon perturbation Eq in Legendre Expansion
	cmb.lmaxnr = 25;		    //Number of relativistic neutrino perturbation Eq in Legendre Expansion
	cmb.nqmax = 15;				//Number of q's at which massive neutrino equations will be written
	cmb.lmaxnu = 25;			//Number of massive neutrino perturbation Eq in Legendre Expansion of each q

    cmb.PerturbationType = 2;   // For tensor alone

    cmb.setparam();


    //printf("\n\n%d %d %d %d\n", cmb.lmaxg,cmb.lmaxnr, cmb.nqmax, cmb.lmaxnu);

    double tau_present = Others.calcconftime(1.0e-8);
    int NumTensor = Others.NumberOfPerturbationVarsTensor();

    //printf("NumPertubationTensor = %d\n", NumTensor);

    double *y;
    y = (double *)malloc((NumTensor+1)*sizeof(double));  // We add +1 because we use y[1:NumScalar]. y[0] is ignored


    FILE *fp;
    fp = fopen("Data_finitialt.d","w");
    for(int l=1;l<=1024;l++)
    {
        cmb.k_mode = l/tau_present;
        cmb.setparam();
        Tensor.finitialt(tau_present, y);
        for(int i=1;i<=NumTensor;i++)
            fprintf(fp,"%e  ",y[i]);
        fprintf(fp, "\n");    
    }
    fclose(fp);
    printf("Data_finitial.d created\n");

    Others.PrintPertubationVarsTensor(y);
}