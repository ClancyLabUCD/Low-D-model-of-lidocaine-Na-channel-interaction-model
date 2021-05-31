//
//  I_Na_TT.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

#ifndef I_Na_TT_h
#define I_Na_TT_h

#include <math.h>
#include "Global_variables.h"
#include "Numerical_Methods.h"

// sodium current from Ten Tusscher model
void Calculate_I_Na_TT(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2004 paper
    
    double a_m, b_m, a_h, b_h, a_j, b_j;
    double tau_m, tau_h, tau_j, m_inf, h_inf, j_inf;
    const double G_Na=14.838; //pS/pF
    
    double E_Na=(R*T/F)*log(Na_out/Cell_ptr->Na_in);
    
    a_m = 1./(1.0 + exp(-60.0-Cell_ptr->V)/5.0);
    b_m = 0.1/(1.0+exp((Cell_ptr->V+35.0)/5.0))+0.1/(1.0+exp((Cell_ptr->V-50.0)/200.0));
    tau_m = a_m*b_m;
    m_inf = 1.0/pow((1.0+exp((-56.86-Cell_ptr->V)/9.03)),2.0);
    
    if(Cell_ptr->V>=-40.0){
        a_h=0.0;
        b_h=0.77/(0.13*(1.0+exp((Cell_ptr->V+10.66)/-11.1)));
        
        a_j=0.0;
        b_j=0.6*exp(0.057*Cell_ptr->V)/(1.0+exp(-0.1*(Cell_ptr->V+32.0)));
    }
    
    else {
        a_h=0.057*exp((80.0+Cell_ptr->V)/-6.8);
        b_h=2.7*exp(0.079*Cell_ptr->V)+3.1E5*exp(0.3485*Cell_ptr->V);
        
        a_j=((-2.5428E4)*exp(0.2444*Cell_ptr->V)-(6.978E-6)*exp(-0.04391*Cell_ptr->V))*(Cell_ptr->V+37.78)/(1+exp(0.311*(Cell_ptr->V+79.23)));
        b_j=0.02424*exp(-0.01052*Cell_ptr->V)/(1.0+exp(-0.1378*(Cell_ptr->V+40.14)));
    }
    
    tau_h = 1/(a_h+b_h);
    tau_j = 1/(a_j+b_j);
    
    h_inf = 1.0/pow((1.0+exp((Cell_ptr->V+71.55)/7.43)),2.0);
    j_inf = h_inf;
    
    //from HH model if time=0 then m0 (resting states) is neglected
    
    Cell_ptr->m = m_inf-(m_inf-Cell_ptr->m)*exp(-dt/tau_m);
    
    Cell_ptr->h = h_inf-(h_inf-Cell_ptr->h)*exp(-dt/tau_h);
    
    Cell_ptr->j = j_inf-(j_inf-Cell_ptr->j)*exp(-dt/tau_j);
    
    Cell_ptr->I_Na = G_Na*(Cell_ptr->V-E_Na)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->m)*(Cell_ptr->h)*(Cell_ptr->j);
}

#endif /* I_Na_TT_h */
