//
//  Ion_conc_dynamics.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

////3-31-18 This code was originally copied from 1D_modified_Ten_Tusscher_with_drug in PLoS_Comp_Bio_2017_paper_C++_codes.

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Functions.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.
////The Functions.h file was split into Non_Na_currents.h, Ion_conc_dynamics.h, and Full_cell_functions.h.

#ifndef Ion_conc_dynamics_h
#define Ion_conc_dynamics_h

#include <math.h>
#include "Global_variables.h"

//New K concentration in myoplasm
void Calculate_K_in(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2004 paper
    double dK_in;
    //I_K_ion_total calculated in Calculate_I_total
    
    dK_in =  -dt*((Cell_ptr->I_K_ion_total + Cell_ptr->I_stim - Cell_ptr->I_axial)*CAP)/(V_cyto*z_K*F);
    Cell_ptr->K_in = Cell_ptr->K_in + dK_in;
}

//New Na concentration in myoplasm
void Calculate_Na_in(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2004 paper
    double dNa_in;
    dNa_in = -1.*dt*(Cell_ptr->I_Na_ion_total*CAP)/(V_cyto*z_Na*F);
    Cell_ptr->Na_in = Cell_ptr->Na_in + dNa_in;
}

//New Ca concentration in cytoplasm
void Calculate_Ca_in(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2006 paper and code.  Form seems to be in Zeng-Rudy 1995.
    double b_Ca, c_Ca, dCa_in;
    const double buffer = 0.2;                //Total cytoplasmic buffer concentration (mM)
    const double K_buffer = 0.001;            // Ca_in half-saturation constant for cytoplasmic buffer (mM)
    
    Cell_ptr->Ca_in_buffer = (Cell_ptr->Ca_in* buffer)/(Cell_ptr->Ca_in + K_buffer);
    dCa_in = dt*((-(Cell_ptr->I_Ca_b+ Cell_ptr->I_p_Ca - 2*Cell_ptr->I_Na_Ca)*CAP/(V_cyto*z_Ca*F))+(V_sr/V_cyto)*(Cell_ptr->I_leak- Cell_ptr->I_up)+ Cell_ptr->I_tr);
    
    b_Ca = buffer - Cell_ptr->Ca_in_buffer - dCa_in - Cell_ptr->Ca_in + K_buffer;
    c_Ca = K_buffer*(Cell_ptr->Ca_in_buffer + dCa_in + Cell_ptr->Ca_in);
    
    Cell_ptr->Ca_in = (sqrt(b_Ca*b_Ca + 4*c_Ca) - b_Ca)/2;
}

//New Ca concentration in SR
void Calculate_Ca_sr(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2006 paper and code.  Form seems to be in Zeng-Rudy 1995.
    double b_jsr, c_jsr, dCa_sr;
    const double buffer_sr = 10;                //Total SR buffer concentration (mM)
    const double K_buffer_sr = 0.3;                // Ca_sr half-saturation constant for SR buffer (mM)
    
    Cell_ptr->Ca_sr_buffer = (Cell_ptr->Ca_sr* buffer_sr)/(Cell_ptr->Ca_sr + K_buffer_sr);
    dCa_sr = dt*(Cell_ptr->I_up - Cell_ptr->I_rel - Cell_ptr->I_leak);
    
    b_jsr = buffer_sr - Cell_ptr->Ca_sr_buffer - dCa_sr - Cell_ptr->Ca_sr + K_buffer_sr;
    c_jsr = K_buffer_sr*(Cell_ptr->Ca_sr_buffer + dCa_sr + Cell_ptr->Ca_sr);
    
    Cell_ptr->Ca_sr = (sqrt(b_jsr*b_jsr + 4*c_jsr) - b_jsr)/2;
}

//New Ca concentration in SS
void Calculate_Ca_ss(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2006 paper and code.  Form seems to be in Zeng-Rudy 1995.
    double b_Ca_ss, c_Ca_ss, dCa_ss;
    const double buffer_ss = 0.4;                //Total SS buffer concentration (mM)
    const double K_buffer_ss = 0.00025;            // Ca_ss half-saturation constant for SR buffer (mM)
    
    Cell_ptr->Ca_ss_buffer = (Cell_ptr->Ca_ss*buffer_ss)/(Cell_ptr->Ca_ss + K_buffer_ss);
    dCa_ss = dt*((-Cell_ptr->I_Ca_L*CAP)/(V_ss*z_Ca*F)+(V_sr/V_ss)* Cell_ptr->I_rel-(V_cyto/V_ss)* Cell_ptr->I_tr);
    
    b_Ca_ss = buffer_ss - Cell_ptr->Ca_ss_buffer - dCa_ss - Cell_ptr->Ca_ss + K_buffer_ss;
    c_Ca_ss = K_buffer_ss*(Cell_ptr->Ca_ss_buffer + dCa_ss + Cell_ptr->Ca_ss);
    
    Cell_ptr->Ca_ss = (sqrt(b_Ca_ss*b_Ca_ss + 4*c_Ca_ss) - b_Ca_ss)/2;
}

#endif /* Ion_conc_dynamics_h */
