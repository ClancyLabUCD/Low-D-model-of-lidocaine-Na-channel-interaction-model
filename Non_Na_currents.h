//
//  Non_Na_currents.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken. All rights reserved.
//

////4-6-18 This code was originally copied from 1D_modified_Ten_Tusscher_with_drug in PLoS_Comp_Bio_2017_paper_C++_codes.

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."  Code in this file was copied from the Functions.h file of the Single_Cell folder of the Moreno source code.  Changes were then made to the code.
////The Functions.h file was split into Non_Na_currents.h, Ion_conc_dynamics.h, and Full_cell_functions.h.

#ifndef Non_Na_currents_h
#define Non_Na_currents_h

#include <math.h>
#include "Global_variables.h"

void Calculate_I_Na_Ca (Cell_param *Cell_ptr, double t){ //From 2006 Ten Tusscher et al. paper
    double k_Na_Ca = 1000.0;            // pA/pF
    /*double k_Na_Ca;
     if (Heart_Failure ==1){k_Na_Ca = 1.65 * 1000;}
     else {k_Na_Ca = 1000;}*/        // used if includeding Heart failure model
    
    double Km_Na = 87.5;            // Na_in half-saturation concentration of NaCa exhanger (mM)
    double Km_Ca = 1.38;            // Ca_in half-saturation concentration of NaCa exhanger (mM)
    double k_Na_Ca_sat = 0.1;        // Saturation factor
    double gamma = 0.35;            // Position of energy barrier controlling voltage dependance of inaca
    double alpha = 2.5;                // Factor enhancing outward nature of I_Na_Ca
    
    Cell_ptr->I_Na_Ca =  k_Na_Ca*(1./(Km_Na*Km_Na*Km_Na + Na_out* Na_out* Na_out))*(1./(Km_Ca+ Ca_out))*
    (1./(1.+k_Na_Ca_sat*exp((gamma-1)* Cell_ptr->V *F/(R*T))))*
    (exp(gamma* Cell_ptr->V *F/(R*T))* Cell_ptr->Na_in* Cell_ptr->Na_in* Cell_ptr->Na_in* Ca_out -
     exp((gamma-1)* Cell_ptr->V *F/(R*T))* Na_out* Na_out* Na_out* Cell_ptr->Ca_in *alpha);
    
}

void Calculate_I_Ca_L(Cell_param *Cell_ptr, double t ){//From the 2006 Ten Tusscher et al. paper
    double d_inf, a_d, b_d, c_d, tau_d;
    double f_inf, a_f, b_f, c_f, tau_f;
    double f2_inf, a_f2, b_f2, c_f2, tau_f2;
    double f_Ca_inf, tau_f_Ca;
    const double G_Ca_L = 3.980E-5;
    /*double G_Ca_L;
     if (Heart_Failure ==1){G_Ca_L = 2.5 * 3.980E-5;}
     else{G_Ca_L = 3.980E-5;}*/ //Used for including heart failure model
    
    d_inf = 1./(1.+exp((-8-Cell_ptr->V)/7.5));
    a_d=1.4/(1.+exp((-35-Cell_ptr->V)/13.0))+0.25;
    b_d=1.4/(1.+exp((Cell_ptr->V+5.)/5.));
    c_d=1./(1.+exp((50.-Cell_ptr->V)/20.));
    
    tau_d = a_d*b_d + c_d;
    Cell_ptr->d = d_inf - (d_inf - Cell_ptr->d)*exp(-dt/tau_d);
    
    f_inf = 1./(1.+exp((Cell_ptr->V+20.)/7.));
    a_f=1102.5*exp(-(Cell_ptr->V+27.)*(Cell_ptr->V+27.)/225.);
    b_f=200./(1.+exp((13.-Cell_ptr->V)/10.));
    c_f=(180./(1.+exp((Cell_ptr->V+30.)/10.)))+20.;
    
    tau_f = a_f + b_f + c_f;
    Cell_ptr->f = f_inf - (f_inf - Cell_ptr->f)*exp(-dt/tau_f);
    
    f2_inf = 0.67/(1.+exp((Cell_ptr->V+35.)/7.))+0.33;
    a_f2= 600.*exp(-(Cell_ptr->V+25.)*(Cell_ptr->V+25.)/170.);
    b_f2= 31./(1.+exp((25.-Cell_ptr->V)/10.));
    c_f2= 16./(1.+exp((Cell_ptr->V+30.)/10.));
    
    tau_f2 = a_f2 + b_f2 + c_f2;
    Cell_ptr->f2 = f2_inf - (f2_inf - Cell_ptr->f2)*exp(-dt/tau_f2);
    
    f_Ca_inf = 0.6/(1.+(Cell_ptr->Ca_ss/0.05)*(Cell_ptr->Ca_ss/0.05))+0.4;
    tau_f_Ca = 80./(1.+(Cell_ptr->Ca_ss/0.05)*(Cell_ptr->Ca_ss/0.05))+2.;
    
    Cell_ptr->f_Ca = f_Ca_inf - (f_Ca_inf - Cell_ptr->f_Ca)*exp(-dt/tau_f_Ca);
    
    Cell_ptr->I_Ca_L = G_Ca_L*Cell_ptr->d* Cell_ptr->f* Cell_ptr->f2* Cell_ptr->f_Ca*4.*(Cell_ptr->V-15.)*(F*F/(R*T))*
    (0.25*exp(2.*(Cell_ptr->V-15.)*F/(R*T))* Cell_ptr->Ca_ss- Ca_out)/(exp(2.*(Cell_ptr->V-15.)*F/(R*T))-1.);
    
}

void Calculate_I_Kr(Cell_param *Cell_ptr, double t) {//From Ten Tusscher et al. 2004 paper
    double a_xr1, b_xr1, xr1_ss, tau_xr1;
    double a_xr2, b_xr2, xr2_ss, tau_xr2;
    const double G_Kr = 0.153;        // nS/pF  in 2006 paper
    double E_Kr = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in);
    
    a_xr1=450./(1.+exp((-45.-Cell_ptr->V)/10.));
    b_xr1=6./(1.+exp((Cell_ptr->V-(-30.))/11.5));
    xr1_ss = 1./(1.+exp((-26.-Cell_ptr->V)/7.));
    tau_xr1 = a_xr1*b_xr1;
    Cell_ptr->xr1 = xr1_ss - (xr1_ss - Cell_ptr->xr1)*exp(-dt/tau_xr1);
    
    a_xr2=3./(1.+exp((-60.-Cell_ptr->V)/20.));
    b_xr2=1.12/(1.+exp((Cell_ptr->V-60.)/20.));
    xr2_ss = 1./(1.+exp((Cell_ptr->V-(-88.))/24.));
    tau_xr2 = a_xr2*b_xr2;
    Cell_ptr->xr2 = xr2_ss - (xr2_ss - Cell_ptr->xr2)*exp(-dt/tau_xr2);
    
    Cell_ptr->I_Kr = G_Kr*sqrt(K_out/5.4)*Cell_ptr->xr1*Cell_ptr->xr2*(Cell_ptr->V-E_Kr);
}

void Calculate_I_Ks(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2006 paper
    double xs_ss, ax_s, bx_s, tau_xs, G_Ks;
    const double PR_NaK = 0.03;
    
    if (Cell_ptr->Cell_type == 1)            // Cell type = (1) for endo; (2) for M,  (3) epi
    {G_Ks = 0.392;}            // nS/pF
    else if (Cell_ptr->Cell_type == 2)
    {G_Ks = 0.098;}            // nS/pF
    else if (Cell_ptr->Cell_type == 3)
    {G_Ks = 0.392;}            // nS/pF        //0.392
    
    double E_Ks = ((R*T)/(z_K*F))*(log((K_out+PR_NaK*Na_out)/(Cell_ptr->K_in + PR_NaK*Cell_ptr->Na_in)));
    
    xs_ss = 1./(1.+exp((-5.-Cell_ptr->V)/14.));
    ax_s =(1400./(sqrt(1.+exp((5.-Cell_ptr->V)/6.))));
    bx_s =(1./(1.+exp((Cell_ptr->V-35.)/15.)));
    tau_xs =(ax_s*bx_s) + 80.;
    
    Cell_ptr->xs = xs_ss - (xs_ss - Cell_ptr->xs)*exp(-dt/tau_xs);
    
    Cell_ptr->I_Ks = G_Ks*Cell_ptr->xs*Cell_ptr->xs*(Cell_ptr->V-E_Ks);
}

void Calculate_I_K1(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2004 paper
    double a_K1, b_K1, K1_s;
    const double G_K1 = 5.405;            // nS/pF
    /*double G_K1;
     if (Heart_Failure ==1){G_K1 = 0.75 * 5.405;}
     else {G_K1 = 5.405;}*/                  //Used if heart failure model is used
    
    double E_K1 = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in);
    
    a_K1 =0.1/(1.+exp(0.06*(Cell_ptr->V-E_K1-200.)));
    b_K1 = (3.*exp(0.0002*(Cell_ptr->V-E_K1+100.)) + exp(0.1*(Cell_ptr->V-E_K1-10.)))/(1.+exp(-0.5*(Cell_ptr->V-E_K1)));
    
    K1_s = a_K1/(a_K1 + b_K1);
    
    Cell_ptr->I_K1=G_K1*K1_s*(Cell_ptr->V-E_K1);
}

void Calculate_I_to(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2004 paper
    double s_inf, tau_s, r_inf, tau_r, G_to;
    if (Cell_ptr->Cell_type == 1)            // Cell type = (1) for endocardial
    {    G_to = 0.073;            // nS/pF
        /*if (Heart_Failure ==1){G_to = 0.64 * 0.073;}
         else {G_to = 0.073;}*/          //Used if heart failure model included
        s_inf = 1./(1.+exp((Cell_ptr->V+28.)/5.));
        tau_s = 1000.*exp(-(Cell_ptr->V+67.)*(Cell_ptr->V+67.)/1000.)+8.;
    }
    
    else if (Cell_ptr->Cell_type == 2)    // Cell type = (2) for M Cell
    {    G_to = 0.294;            // nS/pF
        /*if (Heart_Failure ==1){G_to = 0.64 * 0.294;}
         else {G_to = 0.294;}*/         //Used if heart failure model included
        s_inf = 1./(1.+exp((Cell_ptr->V+20.)/5.));
        tau_s = 85.*exp(-(Cell_ptr->V+45.)*(Cell_ptr->V+45.)/320.)+5./(1.+exp((Cell_ptr->V-20.)/5.))+3.;
    }
    
    else if (Cell_ptr->Cell_type == 3)    // Cell type = (3) for epicardial
    {    G_to = 0.294;            // nS/pF
        /*if (Heart_Failure ==1){G_to = 0.64 * 0.294;}
         else {G_to = 0.294;}*/          //Used if heart failure model included
        s_inf=1./(1.+exp((Cell_ptr->V+20.)/5.));
        tau_s = 85.*exp(-(Cell_ptr->V+45.)*(Cell_ptr->V+45.)/320.)+5./(1.+exp((Cell_ptr->V-20.)/5.))+3.;
    }
    
    double E_to = ((R*T)/(z_K*F))*log(K_out/Cell_ptr->K_in);
    
    r_inf = 1./(1.+exp((20.-Cell_ptr->V)/6.));
    tau_r = 9.5*exp(-(Cell_ptr->V+40.)*(Cell_ptr->V+40.)/1800.)+0.8;
    
    Cell_ptr->r = r_inf-(r_inf-Cell_ptr->r)*exp(-dt/tau_r);
    Cell_ptr->s = s_inf-(s_inf-Cell_ptr->s)*exp(-dt/tau_s);
    
    Cell_ptr->I_to = G_to*Cell_ptr->r*Cell_ptr->s*(Cell_ptr->V-E_to);
}


void Calculate_I_p_Ca(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2004 paper with values from the 2006 paper
    const double G_p_Ca = 0.1238;        // nS/pF
    const double Km_p_Ca = 0.0005;        // Half saturation constant (mM)
    
    Cell_ptr->I_p_Ca=G_p_Ca*(Cell_ptr->Ca_in/(Km_p_Ca + Cell_ptr->Ca_in));
}

void Calculate_I_Kp(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2004 paper
    const double G_Kp = 0.0146;            // nS/pF
    
    double E_Kp=(R*T/F)*log(K_out/Cell_ptr->K_in);
    
    Cell_ptr->I_Kp = G_Kp*(1./(1.+exp((25.-Cell_ptr->V)/5.98)))*(Cell_ptr->V-E_Kp);
}

void Calculate_I_Ca_b(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2004 paper
    const double G_Ca_b =  0.000592;
    
    double E_Ca_b = ((R*T)/(z_Ca*F))*log(Ca_out/Cell_ptr->Ca_in);
    Cell_ptr->I_Ca_b = G_Ca_b*(Cell_ptr->V-E_Ca_b);
}

void Calculate_I_Na_b(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2004 paper
    const double G_Na_b = 0.000290;        // nS/pF
    double E_Na_b=((R*T)/(z_Na*F))*log(Na_out/Cell_ptr->Na_in);
    Cell_ptr->I_Na_b = G_Na_b*(Cell_ptr->V-E_Na_b);
}

void Calculate_I_Na_K(Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2004 paper
    const double I_Na_K_bar = 2.724;    // Maximal I_Na_K (pA/pF)  //From 2006 paper
    
    /*double I_Na_K_bar;
     if (Heart_Failure ==1){I_Na_K_bar = 1.0 * 2.724;}        //0.58 * 2.724
     else {I_Na_K_bar = 2.724;}*/        // Used if heart failure model is included
    
    const double Km_Na_in = 40.;            // Na_in half saturation constant
    const double Km_K_out = 1.;            // K_out half saturation constant
    
    Cell_ptr->I_Na_K= I_Na_K_bar*(K_out*Cell_ptr->Na_in/((K_out+Km_K_out)*(Cell_ptr->Na_in+Km_Na_in)))*
    (1./(1.+0.1245*exp(-0.1*Cell_ptr->V*F/(R*T))+0.0353*exp(-Cell_ptr->V*F/(R*T))));
}


void Calculate_I_rel (Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2006 paper
    double k_Ca_sr, k1_rel, k2_rel, tau_R_bar, R_bar_inf;//dR_bar; //dR_bar is for Ten Tusscher code formulation of R_bar calculation
    const double G_rel = 0.102;                // Max. rate constant of Ca release from JSR due to overload (mM/ms)
    
    const double k1_rel_ = 0.15;            // R to O and RI to I I_rel transition rate (mM^2/ms)
    const double k2_rel_ = 0.045;            // O to I and R to RI I_rel transition rate (mM^2/ms)
    const double k3_rel = 0.060;            // O to R and I to RI I_rel transition rate (ms^-1)
    const double k4_rel = 0.005;            // I to O and RI to I I_rel transition rate (ms^-1)
    const double EC_sr = 1.5;                // Ca_sr half-saturation constant of k_Ca_sr
    const double max_sr = 2.5;                // Max value of k_Ca_sr (dimensionless)
    const double min_sr = 1.;                // Min value of k_Ca_sr (dimensionless)
    
    k_Ca_sr = max_sr-((max_sr-min_sr)/(1.+(EC_sr/Cell_ptr->Ca_sr)*(EC_sr/Cell_ptr->Ca_sr)));
    
    k1_rel = k1_rel_ / k_Ca_sr;
    k2_rel = k2_rel_ * k_Ca_sr;
    
    /*dR_bar = (-k2_rel*Cell_ptr->Ca_ss* Cell_ptr->R_bar + k4_rel*(1.- Cell_ptr->R_bar))*dt;
     Cell_ptr->R_bar = Cell_ptr->R_bar + dR_bar;*/ //This uses Forward Euler and is how it was formulated in Ten Tusscher's code
    
    R_bar_inf = k4_rel/(k2_rel*Cell_ptr->Ca_ss + k4_rel);
    tau_R_bar = 1./(k2_rel*Cell_ptr->Ca_ss + k4_rel);
    
    Cell_ptr->R_bar = R_bar_inf-(R_bar_inf-Cell_ptr->R_bar)*exp(-dt/tau_R_bar); //Analytic solution
    
    Cell_ptr->OO = (k1_rel*Cell_ptr->Ca_ss* Cell_ptr->Ca_ss* Cell_ptr->R_bar)/(k3_rel + k1_rel* Cell_ptr->Ca_ss* Cell_ptr->Ca_ss);
    
    Cell_ptr->I_rel = G_rel*Cell_ptr->OO*(Cell_ptr->Ca_sr - Cell_ptr->Ca_ss);
}

void Calculate_I_up(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2006 paper
    const double G_up = 0.006375;        // mM/ms
    /*double G_up;
     if (Heart_Failure ==1){G_up = 0.64 * 0.006375;}        // 0.64 * 0.006375
     else {G_up = 0.006375;}*/               //Used if a heart failure model is included
    
    const double Km_up = 0.00025;            //mM
    
    Cell_ptr->I_up = G_up*((Cell_ptr->Ca_in* Cell_ptr->Ca_in)/((Cell_ptr->Ca_in* Cell_ptr->Ca_in)+(Km_up*Km_up)));
}

void Calculate_I_leak(Cell_param *Cell_ptr, double t ){//From Ten Tusscher et al. 2006 paper
    const double G_leak = 0.00036;        // mM/ms
    
    /*double G_leak;
     if (Heart_Failure ==1) {G_leak = 0.00025;}    //0.00025
     else {G_leak = 0.00036;}*/      //Used if a heart failure model is included
    
    Cell_ptr->I_leak = G_leak*(Cell_ptr->Ca_sr- Cell_ptr->Ca_in);
}

void Calculate_I_tr( Cell_param *Cell_ptr, double t){//From Ten Tusscher et al. 2006 paper
    const double G_tr = 0.0038;        // mM/ms
    
    Cell_ptr->I_tr = G_tr*(Cell_ptr->Ca_ss- Cell_ptr->Ca_in);
}

#endif /* Non_Na_currents_h */
