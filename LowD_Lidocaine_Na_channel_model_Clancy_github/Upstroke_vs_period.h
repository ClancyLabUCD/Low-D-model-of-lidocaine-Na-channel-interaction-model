//
//  Upstroke_vs_period.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/9/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in "Rate-dependent effects of lidocaine on cardiac dynamics: Development and analysis of a low-dimensional drug-channel interaction model" PLOS Computational Biology
//

//// Copied and edited from PLoS_Comp_Bio_2017_paper_C++_codes/Ten_Tusscher_model_with_drug/Upstroke_vs_period_all_models_mult_time_const.h

#ifndef Upstroke_vs_period_h
#define Upstroke_vs_period_h

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "Global_variables.h"
#include "Non_Na_currents.h"
#include "Ion_conc_dynamics.h"
#include "Na_currents.h"
#include "Full_cell_functions.h"

//This program brings the cell to steady state and then paces it 500 times
void Upstroke_vs_period(){
    Cell_param Cell; //initializing the structure that will hold model variables
    int S1_cycle; //initializing counter for number of stimuli
    
    double t, tTemp, interval; //initializing variable to hold total simulated time, simulated time within a for loop, and the BCL value
    
    std::string dt_val;
    
    if (dt == 0.01) {
        dt_val = "1eneg2";
    }else if (dt == 0.005){
        dt_val = "5eneg3";
    }//setting value of dt to be used in saving the results files
    
    int counter =0; //counts the number of time steps that have been taken
    
    
    std::string Drug_conc; //This will store the drug concentration
    std::string Na_model_str, k_on_fac_str, k_off_fac_str;
    
    double k_on_fac, k_off_fac;
    
    std::ofstream result_file;
    
    for(Na_model = 0; Na_model < 3; Na_model++){//Looping through all Na_models
        switch (Na_model) {
            case 1:{
                Na_model_str = "Simple";
                break;
            }
            case 2:{
                Na_model_str = "Moreno_2011";
                break;
            }
                
            default:{
                Na_model_str = "Ten_Tusscher";
                break;
            }
        }
        t = 0.0;
        tTemp = 0.0; //initializing t and tTemp
        
        //Initilizations for file names
        
        
        counter = 0;
        
        for (int Drug_conc_counter = 0; Drug_conc_counter < 3; Drug_conc_counter++) { //Loops through the various drug concentrations
            
            switch (Drug_conc_counter) {
                case 0:{
                    Drug_conc = "0uM";
                    Drug = 0.0;
                    break;
                }
                case 1:{
                    Drug_conc = "5uM";
                    Drug = 5.0e-6;
                    break;
                }
                case 2:{
                    Drug_conc = "20uM";
                    Drug = 20.0e-6;
                    break;
                }
            }//recording the diffusion value for the file name
            
            if (Na_model == 0) {
                Drug_conc_counter = 3;
            }
            
            k_on_fac = 1.0;
            k_off_fac = 1.0;// possible to alter the drug binding rates with these parameters
            k_on_fac_str = "1";
            k_off_fac_str = "1";
            
                result_file.open((Na_model_str + "_Upstroke_v_Period_dt" + dt_val + "_Lido_conc_" + Drug_conc + "_konfac" + k_on_fac_str + "_kofffac" + k_off_fac_str +  "_stimdur1.txt").c_str());
            
            
                Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
            
            
                for (double jj = 270.0; jj <= 1000.0; jj = jj + 10.0) {//loops through various BCL
                
                    //Loading in initial conditions for each cell in the array
                    Cell.Cell_type = 3; // setting cell type to Epi
                    Cell.V=-86.2;
                    Cell.V_new=-86.2;
                    Cell.dV=0;
                    Cell.Na_in = 7.67;
                    Cell.K_in = 138.3;
                    Cell.Ca_in = 0.00007;
                    Cell.Ca_sr = 1.3;
                    Cell.Ca_ss = 0.00007;
                    Cell.Ca_in_buffer = 0;
                    Cell.Ca_ss_buffer = 0;
                    Cell.Ca_sr_buffer = 0;
                    Cell.m = 0.0;
                    Cell.h = 0.75;
                    Cell.j = 1.0;
                    Cell.b = 0.0;
                
                    //Drug Free States
                    Cell.IC3 = 0;
                    Cell.IC2 = 0;
                    Cell.IF = 0;
                    Cell.IM1 = 0;
                    Cell.IM2 = 0;
                    Cell.C3 = 1;
                    Cell.C2 = 0;
                    Cell.O = 0;
                    Cell.OS =0;
                    //Charged Drug States
                    Cell.DIC3 = 0;
                    Cell.DIC2 =0;
                    Cell.DIF = 0;
                    Cell.DIM1 = 0;
                    Cell.DIM2 = 0;
                    Cell.DC3 = 0;
                    Cell.DC2 =0;
                    Cell.DO = 0;
                    Cell.DOS = 0;
                    Cell.DC1 = 0;
                    //Neutral Drug States
                    Cell.D_IC3 = 0;
                    Cell.D_IC2 = 0;
                    Cell.D_IF = 0;
                    Cell.D_IM1 = 0;
                    Cell.D_IM2 = 0;
                    Cell.D_C3 = 0;
                    Cell.D_C2 = 0;
                    Cell.D_O = 0;
                    Cell.D_OS = 0;
                    Cell.D_C1 = 0;
                
                
                    Cell.C1 = 1 - ( Cell.O + Cell.OS + Cell.C3 + Cell.C2 + Cell.IC3 + Cell.IC2 + Cell.IF + Cell.IM1 + Cell.IM2 + Cell.DO + Cell.DOS + Cell.DC1 + Cell.DC2 + Cell.DC3 + Cell.DIC3 + Cell.DIC2 + Cell.DIF + Cell.DIM1 + Cell.DIM2 + Cell.D_O + Cell.D_OS + Cell.D_C1 + Cell.D_C2 + Cell.D_C3 + Cell.D_IC3 + Cell.D_IC2 + Cell.D_IF + Cell.D_IM1 + Cell.D_IM2);
                
                    /*Cell.mL =  0.00111859;
                     Cell.hL =  0.339310414;*/
                
                    Cell.d = 0 ;
                    Cell.f = 1;
                    Cell.f2 = 1;
                    Cell.f_Ca = 1;
                    Cell.r = 0;
                    Cell.s = 1;
                    Cell.xr1 = 0;
                    Cell.xr2 = 1;
                    Cell.xs = 0;
                    Cell.OO = 0;
                    Cell.R_bar = 1;
                
                    Cell.I_Na = 0;
                    Cell.I_Na_L = 0;
                    Cell.I_Ca_L = 0;
                    Cell.I_Kr = 0;
                    Cell.I_Ks = 0;
                    Cell.I_K1 = 0;
                    Cell.I_Kp = 0;
                    Cell.I_to = 0;
                    Cell.I_Na_Ca = 0;
                    Cell.I_Na_K = 0;
                    Cell.I_p_Ca = 0;
                    Cell.I_Ca_b = 0;
                    Cell.I_Na_b = 0;
                    Cell.I_stim = 0;
                    Cell.I_tr = 0;
                    Cell.I_leak = 0;
                    Cell.I_up = 0;
                    Cell.I_rel = 0;
                
                    Cell.I_Na_ion_total = 0;
                    Cell.I_K_ion_total = 0;
                    Cell.I_Ca_ion_total = 0;
                    Cell.I_total = 0.0;
                    Cell.I_axial=0.0;
                
                    Cell.peak_slope = 0;
                    Cell.V_min = -88.654973;
                    Cell.t_min = 0;
                    Cell.V_thr = -88.654973;
                    Cell.t_thr = 0;
                    Cell.V_max = -88.654973;
                    Cell.t_max = 0;
                    Cell.V_90 = -88.654973;
                    Cell.t_APD90 = 0;
                    Cell.t_APD90_old = 0;
                    Cell.dV_old = 0;
                    Cell.APD90_flag = 0;
                
                    Cell.cycle_num = 0;
                    Cell.AP_flag = 0;
                
                
                    /*********************************************************************************************************************************
                     Begin Time Loop Here
                     *********************************************************************************************************************************/
                
                    for (S1_cycle = 0; S1_cycle <=500; S1_cycle = S1_cycle + 1){                //Enter the number of beats that you want to simulate
                    
                        if (S1_cycle ==0){interval=waitTime;}//allows the model to equilibrate for the waitTime set in Global_variables.h
                    
                        else {interval=jj;}                                                    //sets the BCL
                    
                        Calculate_Reset (&Cell, t, tTemp);//Resets cell specific parameters every beat such as APD, DI, V90 etc.
                        for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){
                            Calculate_Points (&Cell, t); //Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
                        
                            //Updating current calculations for each cell
                            Calculate_I_Na(&Cell, tTemp, k_on_fac, k_off_fac); //This is where the Na channel model is chosen
                            Calculate_I_Ca_L(&Cell, tTemp );
                            Calculate_I_Kr(&Cell, tTemp );
                            Calculate_I_Ks(&Cell, tTemp );
                            Calculate_I_K1(&Cell, tTemp );
                            Calculate_I_Kp(&Cell, tTemp );
                            Calculate_I_to(&Cell, tTemp );
                            Calculate_I_Na_Ca(&Cell, tTemp );
                            Calculate_I_Na_K(&Cell, tTemp );
                            Calculate_I_p_Ca(&Cell, tTemp );
                            Calculate_I_Ca_b(&Cell, tTemp );
                            Calculate_I_Na_b(&Cell, tTemp );
                            Calculate_I_total(&Cell, t, tTemp);
                        
                            //Updating ionic concentrations
                            Calculate_Na_in(&Cell, tTemp );
                            Calculate_K_in(&Cell, tTemp );
                            Calculate_I_tr(&Cell, tTemp);
                            Calculate_I_leak(&Cell, tTemp);
                            Calculate_I_up(&Cell, tTemp);
                            Calculate_I_rel(&Cell, tTemp);
                            Calculate_Ca_sr(&Cell, tTemp );
                            Calculate_Ca_ss(&Cell, tTemp );
                            Calculate_Ca_in(&Cell, tTemp );
                        
                            //Calculates t_min, V_min, t_max, V_max etc. etc.; Update calculates axial currents and updates the voltage
                            Calculate_Update (&Cell, t, tTemp);//Calculate new V
                        
                        
                            t = t+dt; //update time
                        
                            /*******************************************************************************************************************************************************
                             Beginning of result files
                             *******************************************************************************************************************************************************/
                        
                        
                            if (((interval - tTemp)<dt) && S1_cycle > 498) {
                                result_file << jj << ", " << S1_cycle << ", " << Cell.peak_slope << ", " << Cell.APD_90 << ", " << Cell.b_star << std::endl;//recording the BCL, number of stimuli, peak Upstroke, APD90, and fraction of channels bound to drug during the peak upstroke velociey for each stimuli
                            }
                        
                        
                            /*******************************************************************************************************************************************************
                             End of result files
                             *******************************************************************************************************************************************************/
                        
                            counter++;
                        
                            if (fmod(t,5000)<dt) {std::cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval << ", " << "Drug = "<< Drug_conc << std::endl;}//outputting progress to command line
                        
                        }//End of For tTemp ==
                    
                    }//End of for S1
                
                } // End BCL loop
            
                result_file.close();
        } //End Concentration loop
    } // End Drug model loop
}

#endif /* Upstroke_vs_period_h */
