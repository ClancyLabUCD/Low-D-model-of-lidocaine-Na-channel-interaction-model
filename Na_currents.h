//
//  Na_currents.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

#ifndef Na_currents_h
#define Na_currents_h

#include <math.h>
#include "Global_variables.h"
#include "WT_Lido_Implicit.h"
#include "I_Na_Simple_Lidocaine.h"
#include "I_Na_TT.h"


//This function uses the correct I_Na model depending on the value of Na_model
void Calculate_I_Na(Cell_param *Cell_ptr, double t, double k_on_fac, double k_off_fac){
    if (Na_model == 0){
        Calculate_I_Na_TT(Cell_ptr,t);
    } else if (Na_model == 1) {
        Calculate_I_Na_Simple_Lidocaine(Cell_ptr, t, k_on_fac, k_off_fac);
    } else if (Na_model == 2){
        WT_SCN5A_Lidocaine_function(Cell_ptr, t, k_on_fac, k_off_fac); //This is taken from the code for Moreno et al. 2011
    }
}

#endif /* Na_currents_h */
