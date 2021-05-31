//
//  Numerical_Methods.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

////4-6-18 This code was originally copied from 1D_modified_Ten_Tusscher_with_drug in PLoS_Comp_Bio_2017_paper_C++_codes.

#ifndef Numerical_Methods_h
#define Numerical_Methods_h


#include <math.h>
#include "Global_variables.h"

//Rush Larsen for normal Hodgkin-Huxley gates of the Na channel.  This function takes the initial value for the gating variable (x_init), the steady state value (x_inf), and the time constant (tau_x) for the given voltage and use the analytic solution for the gating variable (assuming V is constant) to update the gating variable.
double Rush_Larsen(double x_init, double x_inf, double tau_x){
    double x_final;
    x_final = x_inf-(x_inf-x_init)*exp(-dt/tau_x);
    
    return x_final;
}

/*
//Code for updating a gate immobilization binding scheme.  This function takes the initial value (_init), the stead state value (_inf), and the time constant (tau) for b and the gating variable for the state drug cannot bind to (x).  For example, for inactive state binding drugs, x = h, but for drugs that bind to noninactive state channels, x = hbar.  The finally input for this function is the pointer to the array that the updated state variables will be placed in.
void Gate_Immob(double x_init, double x_inf, double tau_x, double b_init, double b_inf, double tau_b, double *vars_final){
    
    vars_final[0] = (1.0-b_init)*x_inf -((1.0-b_init)*x_inf - x_init)*exp(-dt/tau_x); //This is the line to update the gating variable for the state drug cannot bind to.
    
    vars_final[1] = (1.0-x_init)*b_inf -((1.0-x_init)*b_inf - b_init)*exp(-dt/tau_b); //This is the line to update b (the fraction of channels bound to drug)
}
 */

#endif /* Numerical_Methods_h */
