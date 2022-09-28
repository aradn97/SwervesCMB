#ifndef VPC_H_
#define VPC_H_

#include<stdio.h>
#include<math.h>
#include "variables.h"
#include "initval.h"

// You need to allocate the functions here
class VPC{

private :

/****************************************************************************************/
//   These variables are calculated from the recombination process							 //
/****************************************************************************************/	
double alpha = 0;// 1.8; beta = -1.8
double beta = -alpha;
double testflag = 0;
public :

double myadotoa(double a)
{
    double grho2 = Variables.grhom * (Variables.omegac + Variables.omegab) * a + Variables.grhog + Variables.grhor * Variables.annur + Variables.grhom * Variables.omegavdyn * pow(a, 4) + Variables.grhom*Variables.omegak*a*a; 
    return sqrt(grho2 / 3.0)/a/a; 
}

double Variable_C(double a)
{
    if(Variables.VPCflag == 1)
    {
        if(testflag == 1)
            return 0.379;
        return exp(pow(a,alpha)-1);
        }
    return 1;
    }
double Variable_ionization(double a)
{
    if(Variables.VPCflag == 1)
        return 1;  //exp(0.75*(pow(a,alpha)-1));   // If ionization energy changes. 
    return 1;
    }    
double Variable_G(double a)
{
    if(Variables.VPCflag == 1)   
    {
        if(testflag == 1)
            return 0.0498;     
        return exp(3.0*(pow(a,alpha)-1));
        }
    return 1;    
    }
double Variable_h(double a)
{
    if(Variables.VPCflag == 1)     
    {   
        if(testflag == 1)
            return 0.379;        
        return exp(pow(a,alpha)-1);
        }
    return 1;    
    }
double Variable_mu(double a)
{
    if(Variables.VPCflag == 1)   
    {     
        if(testflag == 1)
            return 0.379;                
        return exp(pow(a,alpha)-1);  
        }
    return 1;        
    }
double Variable_kB(double a)
{
    if(Variables.VPCflag == 1)   
    {   
        if(testflag == 1)
            return 0.286;                  
        return exp(1.25*(pow(a,alpha)-1));  
        }
    return 1;
    }
double Variable_Lambda(double a)
{
    if(Variables.VPCflag == 1)  
    {  
        if(testflag == 1)
            return 1;                          
        return exp(pow(a,alpha)-1);  
        }
    return 1;
    }
double Variable_sigmaTh(double a)
{
    return 1;
    }


double Variable_Cdot_o_C(double a)
{
    if(Variables.VPCflag == 1)    
    {
        if(testflag == 1)
            return 0;             
        return pow(a,alpha)*(alpha) * myadotoa(a);    
        }
    return 0;
    }
double Variable_hdot_o_h(double a)
{
    if(Variables.VPCflag == 1)   
    { 
        if(testflag == 1)
            return 0;                     
        return pow(a,alpha+1)*(alpha-1) + 1;    
        }
    return 0;
    }
double Variable_Gdot_o_G(double a)
{
    if(Variables.VPCflag == 1)        
        return 3*Variable_Cdot_o_C(a);
    return 0;
    }
double Variable_mudot_o_mu(double a)
{
    if(Variables.VPCflag == 1)  
    {  
        if(testflag == 1)
            return 0;                     
        return pow(a,alpha+1)*(alpha-1) + 1.0;
        }
    return 0;
    }
double Variable_kBdot_o_kB(double a)
{
    if(Variables.VPCflag == 1)    
    {   
        if(testflag == 1)
            return 0;                      
        return 1.25*(pow(a,alpha+1)*(alpha-1) + 1.0);
        }
    return 0;
    }
double Variable_Lambdadot_o_Lambda(double a)
{
    if(Variables.VPCflag == 1)    
    {
        if(testflag == 1)
            return 0;                     
        return pow(a,-alpha)*(-alpha) * myadotoa(a);    
        }
    return 0;
    }

double Variable_G_c(double a)
{
    if(Variables.VPCflag == 1)        
        return Variable_Gdot_o_G(a) - 2*Variable_Cdot_o_C(a);
    return 0;
    }
 

double Variable_grhom(double a)
{
    if(Variables.VPCflag == 1)        
        return Variables.grhom * Variable_G(a) / pow(Variable_C(a),2);
    return Variables.grhom;
    }

double Variable_grhog(double a)
{
    if(Variables.VPCflag == 1)        
        return Variables.grhog * Variable_G(a) / pow(Variable_C(a),2);
    return Variables.grhog;
    }

double Variable_grhor(double a)
{
    if(Variables.VPCflag == 1)        
        return Variables.grhor * Variable_G(a) / pow(Variable_C(a),2);
    return Variables.grhor;
    }

double Variable_grhonr(double a)
{
    if(Variables.VPCflag == 1)        
        return Variables.grhonr * Variable_G(a) / pow(Variable_C(a),2);
    return Variables.grhonr;    
    }

double Variable_G_c_rho_Lambda(double a, double rho)
{
    double energy;
    initval Initval;
    energy = Variable_G_c(a)*rho+pow(Variable_C(a),4)/8/Initval.pi/Variable_G(a)*Variable_Lambda(a);
    if(Variables.VPCflag == 1)
        return energy;
    return 0;
    }

};

#endif