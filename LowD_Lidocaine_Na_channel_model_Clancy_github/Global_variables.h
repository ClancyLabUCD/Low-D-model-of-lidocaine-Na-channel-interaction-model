//
//  Global_variables.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in "Rate-dependent effects of lidocaine on cardiac dynamics: Development and analysis of a low-dimensional drug-channel interaction model" PLOS Computational Biology
//

////4-6-18 This code was originally copied from 1D_modified_Ten_Tusscher_with_drug in PLoS_Comp_Bio_2017_paper_C++_codes.

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Global_variables.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.

#ifndef Global_variables_h
#define Global_variables_h


#include <math.h>
/*********************************************************************************************************
 Universal Constants
 *********************************************************************************************************/

//Ion Valences and Universal Constants
const double R = 8314.472;                // mJ/mol*K
const double T = 310;                    // K
const double F = 96485.3415;            // C/mol
const double Cm = 2.0;                    // uF/ cm^2
const double CAP = 0.185;                //Cellular capacitance  From Ten Tusscher code
const double rho = 162.;                    // ohm*cm
const double z_Na = 1.;
const double z_Ca = 2.;
const double z_K = 1.;
const double dx= 0.01;                    //in cm

//Cell Geometry
const long double pi = 3.141592;
const double S_cg = 0.2;                //Surface to volume ratio (um^-1)
const double Diffusion = 0.00154;       //Diffusion coefficient for V from ten Tusscher et al. 2004

//Intracellular volumes
const double V_cyto=0.016404;
const double V_sr=0.001094;
const double V_ss=0.00005468; //Not sure what these units are, since they don't match with what is in the paper, but do match what is in the code.

//Extraceullular concentrations (mM)
const double Na_out = 140;
const double Ca_out = 2.0;
const double K_out = 5.4;

//Drug parameters
double Drug = 0.0; //Concentration of drug (M^-1)

//Beating parameters
const double I_duration = 1.0;
const double stimulus = -80;
const double dt = 0.005; //0.01; //Normally dt set to 0.01, but to check convergence, dt set to .005
const double waitTime = 10000; //10 seconds  //30*60000; //30 minutes for I.C.                //If sim_type =0, waitTime ==30000; if sim_type =1, waitTime==0.01;

int Na_model = 1; //Determines whether TT (0), the Simple model (1), or the Moreno 2011 model (2) is used for I_Na



/*********************************************************************************************************
 Structures
 *********************************************************************************************************/
struct Cell_param{ //Moreno Code uses typedef to define this.
    double V, dV, V_new;
    double Na_in, K_in, Ca_in, Ca_sr, Ca_ss, Ca_in_buffer, Ca_ss_buffer, Ca_sr_buffer;
    double I_Na, I_Na_L, I_Ca_L, I_Kr, I_Ks, I_K1, I_Kp, I_to, I_Na_Ca, I_Na_K, I_p_Ca, I_Ca_b, I_Na_b, I_stim;
    double I_Na_ion_total, I_Ca_ion_total, I_K_ion_total, I_total, I_axial;
    double I_tr, I_leak, I_up, I_rel;
    int Cell_type;
    double m, h, j, b, i_O, d, f, f2, f_Ca, r, s, xr1, xr2, xs, OO, R_bar; //b was added (which is the fraction of channels bound to drug), and i_O (which is the fraction of channels channels in the non-inactivated, non-drug bound state in the Full Guarded Receptor model)
    
    double O , OS , C1 , C2 , C3 , IC3 , IC2 , IF , IM1 , IM2 ;
    double DO , DOS , DC1 , DC2 , DC3 , DIC3 , DIC2 , DIF , DIM1 , DIM2 ;
    double D_O , D_OS , D_C1 , D_C2 , D_C3 , D_IC3 , D_IC2 , D_IF , D_IM1 , D_IM2 ;
    
    double O_n , OS_n , C1_n , C2_n , C3_n , IC3_n , IC2_n , IF_n , IM1_n , IM2_n ;
    double DO_n , DOS_n , DC1_n , DC2_n , DC3_n , DIC3_n , DIC2_n , DIF_n , DIM1_n , DIM2_n ;
    double D_O_n , D_OS_n , D_C1_n , D_C2_n , D_C3_n , D_IC3_n , D_IC2_n , D_IF_n , D_IM1_n , D_IM2_n ;
    
    double O_o , OS_o , C1_o , C2_o , C3_o , IC3_o , IC2_o , IF_o , IM1_o , IM2_o ;
    double DO_o , DOS_o , DC1_o , DC2_o , DC3_o , DIC3_o , DIC2_o , DIF_o , DIM1_o , DIM2_o ;
    double D_O_o , D_OS_o , D_C1_o , D_C2_o , D_C3_o , D_IC3_o , D_IC2_o , D_IF_o , D_IM1_o , D_IM2_o ;
    
    //double mL, hL; //These variables are for a heart failure model not used here
    
    double peak_slope, t_min, V_min, t_thr, V_thr, t_thr_old, t_max, V_max, I_Na_peak, t_I_Na_peak, I_Na_L_peak, t_I_Na_L_peak, b_star;
    double t_APD90, V_90, dV_old;
    int APD90_flag, cycle_num, AP_flag;
    double APD_90, DI;
    double t_APD90_old;
    double CV_1cell, CV_2cell, CV_10cell, T; //CV_1cell and CV_2cell calculate CV based on time it takes AP to travel 1 or 2 cells (to see if there is a difference in the estimation)
    
    double I_total_old, dI_total, dI_total_old;
};


#endif /* Global_variables_h */
