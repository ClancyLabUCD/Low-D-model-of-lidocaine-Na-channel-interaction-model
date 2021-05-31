//
//  I_Na_Simple_Lidocaine.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

#ifndef I_Na_Simple_Lidocaine_h
#define I_Na_Simple_Lidocaine_h

#include <math.h>
#include "Global_variables.h"
#include "Numerical_Methods.h"

////3-31-18 This code was originally copied from 1D_modified_Ten_Tusscher_with_drug in PLoS_Comp_Bio_2017_paper_C++_codes.
//// Note: G_Na was adjusted so that Conduction Velocity matches that of Moreno et al. 2011 with no drug at a BCL of 1000ms.


//My simple model of the sodium current.  These parameters were found on 6-14-16. h parameters come from taking the values from Moreno Methods paper and then fitting h_infty (as a function of 2 parameters) to SS Availability and then holding these parameters constant while tau_h is fitted to tau_50.  Then, holding the h parameters constant, the parameters for m were fit to SS Activation (here the full Na-channel model was used, not just m) and tau_m data. Parameters were then adjusted to 310K
void Calculate_I_Na_Simple_Lidocaine(Cell_param *Cell_ptr, double t, double k_on_fac, double k_off_fac){
    
    double a_m, b_m, a_h, b_h;
    double tau_m, tau_h, m_inf, h_inf;
    const double G_Na=20.0;//adjusted to match CV of Moreno et al. 2011 with no drug at a BCL of 1000ms.  originally 17.25; //nS/pF
    
    double E_Na=(R*T/F)*log(Na_out/Cell_ptr->Na_in);
    
    a_m = 45.43*exp((Cell_ptr->V)/13.78);
    b_m = 0.6628*exp((Cell_ptr->V)/-23.25);
    
    a_h = 6.169e-05*exp((Cell_ptr->V)/-9.328);
    b_h = 14.15*exp((Cell_ptr->V)/14.91); //m params adjusted for temperature.
    
    m_inf = a_m/(a_m + b_m);
    h_inf = a_h/(a_h + b_h);
    
    tau_m = 1.0/(a_m + b_m);
    tau_h = 1.0/(a_h + b_h);
    
    //double k_on = 979; //M^-1ms^-1
    //double K_D_0 = 3.28e-6; //M^-1
    //double k_off = k_on*K_D_0;
    
    const double diffusion=500;
    const double kon_open = diffusion; //on rate for binding to the open state in the Moreno et al. 2011 model.
    double k_on=kon_open/2*k_on_fac;
    double k_off=3.4*(1e-6)*diffusion*k_off_fac; //Moreno et al. 2011 rate constants.
    
    const double pH=7.4;
    const double pKa=7.6;
    const double portion = 1/(1+ pow(10, (pH-pKa)) );
    
    const double drug_neutral=Drug*(1-portion); //Calculating the portion of Liodocaine that is neutral at physiological pH
    
    double tau_b = 1.0/((1-(Cell_ptr->h))*drug_neutral*k_on + k_off);
    double b_inf = (1-(Cell_ptr->h))*drug_neutral*k_on/((1-(Cell_ptr->h))*drug_neutral*k_on + k_off);
    //current values of tau_b and b_inf based on current value of h
    
    Cell_ptr->m = Rush_Larsen(Cell_ptr->m, m_inf, tau_m);
    Cell_ptr->h = Rush_Larsen(Cell_ptr->h, h_inf, tau_h);
    Cell_ptr->b = Rush_Larsen(Cell_ptr->b, b_inf, tau_b);
    
    /*
    //This is for the Full gate immobilization binding scheme with inactive state binding
    double vars_final[2]; //Array that will hold the new h and b values
    Gate_Immob(Cell_ptr->h, h_inf, tau_h, Cell_ptr->b, b_inf, tau_b, vars_final);
    
    Cell_ptr->h = vars_final[0];
    Cell_ptr->b = vars_final[1];//Saving the new h and b values
    */ //Using the simplified form of the model for when drug binding is much slower than inactivation
    
    Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(1 - (Cell_ptr->b));
    
}

#endif /* I_Na_Simple_Lidocaine_h */
