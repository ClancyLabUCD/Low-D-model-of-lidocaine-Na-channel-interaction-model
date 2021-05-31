//
//  Full_cell_functions.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Functions.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.
////The Functions.h file was split into Non_Na_currents.h, Ion_conc_dynamics.h, and Full_cell_functions.h.

#ifndef Full_cell_functions_h
#define Full_cell_functions_h


#include <math.h>
#include "Global_variables.h"

void Calculate_Reset (Cell_param *Cell_ptr, double t, double tTemp){
    //resetting variables that store cell characteristics (like peak upstroke velocity or APD90) and variables used to calculate these characteristics.
    Cell_ptr->t_APD90_old = Cell_ptr->t_APD90;
    Cell_ptr->t_APD90 = 0.0;
    Cell_ptr->peak_slope = 0.0;
    Cell_ptr->b_star = 0.0;
    Cell_ptr->CV_1cell = 0.0;
    Cell_ptr->CV_2cell = 0.0;
    Cell_ptr->CV_10cell = 0.0;
    Cell_ptr->V_min = Cell_ptr->V;
    Cell_ptr->V_max = -40;
    Cell_ptr->t_thr_old = Cell_ptr->t_thr;
    
    Cell_ptr->V_thr = 0.0;
    Cell_ptr->DI = 0.0;
    Cell_ptr->T = 0.0;
    Cell_ptr->AP_flag = 0; //resetting flag for the fact that V went above -40mV
    
    
}

void Calculate_Update (Cell_param *Cell_ptr, double t, double tTemp){
    Cell_ptr->dV_old = Cell_ptr->dV;
    Cell_ptr->dV = -1*(Cell_ptr->I_total)*dt;
    Cell_ptr->V = Cell_ptr->V + Cell_ptr->dV;
}

void Calculate_I_total(Cell_param *Cell_ptr, double t, double tTemp){
    Cell_ptr->I_total_old = Cell_ptr->I_total;
    Cell_ptr->dI_total_old = Cell_ptr->dI_total;
    
    if ((t>waitTime) && (tTemp <=I_duration)) {Cell_ptr->I_stim=stimulus;}
    else {Cell_ptr->I_stim=0;}
    
    Cell_ptr->I_Na_ion_total = Cell_ptr->I_Na + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;
    /*Cell_ptr->I_Na_ion_total = Cell_ptr->I_Na + Cell_ptr->I_Na_L + Cell_ptr->I_Na_b + 3*Cell_ptr->I_Na_K + 3*Cell_ptr->I_Na_Ca;*/ //I_Na_L used if heart failure model included.
    Cell_ptr->I_K_ion_total = Cell_ptr->I_Kr + Cell_ptr->I_Ks + Cell_ptr->I_K1 + Cell_ptr->I_Kp - 2*Cell_ptr->I_Na_K + Cell_ptr->I_to;
    Cell_ptr->I_Ca_ion_total = Cell_ptr->I_Ca_L + Cell_ptr->I_p_Ca + Cell_ptr->I_Ca_b - 2*Cell_ptr->I_Na_Ca;
    Cell_ptr->I_total = Cell_ptr->I_stim + Cell_ptr->I_Na_ion_total + Cell_ptr->I_K_ion_total + Cell_ptr->I_Ca_ion_total;
    
    Cell_ptr->dI_total = (Cell_ptr->I_total - Cell_ptr->I_total_old)/dt;
}

void Calculate_Points(Cell_param *Cell_ptr, double t){
    
    if ((Cell_ptr->AP_flag == 0) && (Cell_ptr->V > -40.0)) {
        Cell_ptr->AP_flag = 1; //updating flag that V went above -40mV
    }
    
    //Threshold voltage point catcher needed for t_min
    if ((Cell_ptr->AP_flag == 1) && (Cell_ptr->dV/dt) > Cell_ptr->peak_slope) {
        
        Cell_ptr->peak_slope = Cell_ptr->dV/dt; //recording peak slope
        Cell_ptr->V_thr = Cell_ptr->V; //recording V of peak slope
        Cell_ptr->t_thr = t; //recording time of peak slope
        
        if(Na_model == 1){ //recording fraction of channels bound to drug during peak upstroke
            Cell_ptr->b_star = Cell_ptr->b;
        } else if(Na_model == 2){
            Cell_ptr->b_star = Cell_ptr->DIC3 + Cell_ptr->DIC2 + Cell_ptr->DIF + Cell_ptr->DIM1 + Cell_ptr->DIM2 + Cell_ptr->DC3 + Cell_ptr->DC2 + Cell_ptr->DO + Cell_ptr->DOS + Cell_ptr->DC1 + Cell_ptr->D_IC3 + Cell_ptr->D_IC2 + Cell_ptr->D_IF + Cell_ptr->D_IM1 + Cell_ptr->D_IM2 + Cell_ptr->D_C3 + Cell_ptr->D_C2 + Cell_ptr->D_O + Cell_ptr->D_OS + Cell_ptr->D_C1;
        }
        
        //Diastolic Interval
        Cell_ptr->DI = Cell_ptr->t_thr - Cell_ptr->t_APD90;
        //Period
        Cell_ptr->T = Cell_ptr->t_thr - Cell_ptr->t_thr_old;
    }
    
    if ((Cell_ptr->dV_old < 0.0)&&(Cell_ptr->dV > 0.0)&&(Cell_ptr->V < Cell_ptr->V_min)) {
        
        Cell_ptr->V_min = Cell_ptr->V; //recording minimum V during diastolic interval
        Cell_ptr->t_min = t; //recording time of minimum V during diastolic interval
        
    }
    
    if ((Cell_ptr->AP_flag == 1) && (Cell_ptr->dV_old > 0.0)&&(Cell_ptr->dV < 0.0)&&(Cell_ptr->V > Cell_ptr->V_max)) {
        
        Cell_ptr->V_max = Cell_ptr->V;  //recording maximum V
        Cell_ptr->t_max = t; //recording time of max V
        
        Cell_ptr->V_90 = Cell_ptr->V_min + 0.1*(Cell_ptr->V_max - Cell_ptr->V_min); //calculating potential at which the cell is 90% repolarized for APD90
        
        Cell_ptr->APD90_flag = 1; //stating APD90 can now be calculated.
    }
    
    if ((Cell_ptr->V < Cell_ptr->V_90)&&(Cell_ptr->APD90_flag == 1)) {
        
        Cell_ptr->t_APD90 = t; //saving time of APD90
        Cell_ptr->APD_90 = t - Cell_ptr->t_thr; //Calculating APD90
        
        Cell_ptr->APD90_flag = 0; //ensures APD90 is not recalculated
        Cell_ptr->cycle_num = Cell_ptr->cycle_num + 1; //increases the counter for the period the cell is on by one.
        
    }
    
}

#endif /* Full_cell_functions_h */
