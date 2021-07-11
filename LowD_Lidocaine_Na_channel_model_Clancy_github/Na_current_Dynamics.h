//
//  Na_current_Dynamics.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in "Rate-dependent effects of lidocaine on cardiac dynamics: Development and analysis of a low-dimensional drug-channel interaction model" PLOS Computational Biology
//

#ifndef Na_current_Dynamics_h
#define Na_current_Dynamics_h


#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
//#include <omp.h>
#include <unistd.h>
#include "Global_variables.h"
#include "Non_Na_currents.h"
#include "Ion_conc_dynamics.h"
#include "Na_currents.h"
#include "Full_cell_functions.h"

//This program brings the cell to steady state and then paces it 500 times
void Na_current_Dynamics(){
    Cell_param Cell;
    int S1_cycle;
    
    double t, tTemp, interval, BCL;
    
    double Drug_bound = 0.0; //This variable will hold the fraction of channels bound to drug
    double k_on_fac, k_off_fac; //initializing the factors for increasing k_on and k_off for sensitivity checks
    
    
    std::string dt_val, BCL_val;
    
    if (dt == 0.01) {
        dt_val = "1eneg2";
    }else if (dt == 0.005){
        dt_val = "5eneg3";
    }//setting value of dt to be used in saving the results files
    
    int counter = 0;
    
    std::ofstream result_file, result_file0, result_file2, result_file3, result_file4;
    std::string Na_model_str, Drug_conc, k_on_fac_str, k_off_fac_str;
    
    for (int BCL_counter = 0; BCL_counter < 3; BCL_counter++) {
        if (BCL_counter == 0) {
            BCL = 300.0;
            BCL_val = "300";
        } else if (BCL_counter == 1) {
            BCL = 750.0;
            BCL_val = "750";
        } else if (BCL_counter == 2) {
            BCL = 1000.0;
            BCL_val = "1000";
        }
    
    for (Na_model = 0; Na_model < 3; Na_model++) {
        
        if (Na_model == 0) {
            Na_model_str = "Ten_Tusscher";
        } else if(Na_model == 1){
            Na_model_str = "Simple";
        } else if(Na_model == 2){
            Na_model_str = "Moreno_2011";
        }
        
        t = 0.0;
        tTemp = 0.0; //initializing t and tTemp
        
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
            
            
                
        
                result_file.open((Na_model_str+ "_Na_current_Dynamics_dt" + dt_val + "_Lido_conc_" + Drug_conc + "_konfac" + k_on_fac_str + "_kofffac" + k_off_fac_str + "_BCL" + BCL_val + "_stimdur1.txt").c_str());            // This file contains all of the voltages for each cell (the master-esque file) for the final 5 beats
                result_file0.open((Na_model_str+ "_Na_current_Dynamics_dt" + dt_val + "_Lido_conc_" + Drug_conc + "_konfac" + k_on_fac_str + "_kofffac" + k_off_fac_str +"_BCL" + BCL_val + "_stimdur1_param.txt").c_str()   );    // These files hold data for each cell cycle (to compare V_min, V_max, V_90 etc. to compare between cells)
                result_file3.open((Na_model_str+ "_Na_current_Dynamics_dt" + dt_val + "_Lido_conc_" + Drug_conc + "_konfac" + k_on_fac_str + "_kofffac" + k_off_fac_str +"_BCL" + BCL_val + "_stimdur1_IC_check.txt").c_str()   );    // This file contains all the cell variables at two points to check if steady state was reached.
                result_file4.open((Na_model_str+ "_Na_current_Dynamics_dt" + dt_val + "_Lido_conc_" + Drug_conc + "_konfac" + k_on_fac_str + "_kofffac" + k_off_fac_str +"_BCL" + BCL_val + "_stimdur1_initial_APs.txt").c_str()   );    // This file contains all the variables in the result_file, but for the first 50 APs
    
                Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
        
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
        
    
    
    
                //Cell.Cell_type = 3;        // Cell type = (1) for endo; (2) for M,  (3) epi
                /*********************************************************************************************************************************
                 Begin Time Loop Here
                 *********************************************************************************************************************************/
    
                t = 0.0; //resetting t
            
                for (S1_cycle = 0; S1_cycle <=500; S1_cycle = S1_cycle + 1){                    //Enter the number of beats that you want to simulate
        
                    if (S1_cycle ==0){interval=waitTime;}
        
                    else {interval=BCL;} //This is a frequency of 80 BPM                                                    //Enter the BCL of the simulation here
        
                
                    //Resets cell specific parameters every beat such as APD, DI, V90 etc.
                    Calculate_Reset (&Cell, t, tTemp);
                    for (tTemp = 0; tTemp <= interval; tTemp=tTemp+dt){
                        Calculate_Points (&Cell, t); //Check points before updating to the next time point
            
                        //Updating current calculations for each cell
                        Calculate_I_Na(&Cell, tTemp, k_on_fac, k_off_fac ); //This is where the Na channel model is chosen
                        //Calculate_I_Na_L(&Cell, t, tTemp ); //Used if heartfailure model included
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
                        Calculate_Update (&Cell, t, tTemp);
                
                
                        t = t+dt; //update time
                
                        /*if (sim_type ==0){                                        //This will rewrite state file name for initial conditions (a very long hold time: 30,000ms)
                         if (S1_cycle ==0 && tTemp >= (interval-dt)  ){
                         FILE *fp = fopen (stateFileName, "w" );
                         fwrite (&Cell, sizeof(Cell), 1, fp);
                         fclose (fp);
                         exit(0);}
                         }
             
                         if (sim_type ==1){                                        //This will rewrite state file name2 for initial conditions with drug
                         if (S1_cycle ==5000 && tTemp >= (interval-dt)  ){
                         FILE *fp = fopen (stateFileName2, "w" );
                         fwrite (&Cell, sizeof(Cell), 1, fp);
                         fclose (fp);
                         exit(0);}
                         }*/
                
            
            
            
                    /*******************************************************************************************************************************************************
                    Beginning of result files
                 *******************************************************************************************************************************************************/
            
            
            
                        //This file is for individual cells that calculates all of the min, max, V_90 parameters.
            
                        if (counter%10==0 && S1_cycle >=495 ){
                            if(Na_model == 1){
                                Drug_bound = Cell.b;
                            } else if(Na_model == 2){
                                Drug_bound = Cell.DIC3 + Cell.DIC2 + Cell.DIF + Cell.DIM1 + Cell.DIM2 + Cell.DC3 + Cell.DC2 + Cell.DO + Cell.DOS + Cell.DC1 + Cell.D_IC3 + Cell.D_IC2 + Cell.D_IF + Cell.D_IM1 + Cell.D_IM2 + Cell.D_C3 + Cell.D_C2 + Cell.D_O + Cell.D_OS + Cell.D_C1;
                            }
                            result_file << t << ", " << Cell.V << ", " << Cell.I_Na<< ", " << Cell.I_Ca_L << ", " << Cell.I_Na_L << ", "<< Cell.I_up << ", "<< Cell.I_leak << ", "  << Cell.I_Ks << ", " << Cell.Ca_in << ", " << Cell.Ca_ss << ", " << Drug_bound << std::endl;}
                
                        if ((S1_cycle>0)&&(tTemp >= interval - dt)) {
                            result_file0 << Cell.t_min << ", " << Cell.V_min << ", " << Cell.t_thr<< ", " << Cell.V_thr << ", " << Cell.t_max << ", " << Cell.V_max << ", "<< Cell.t_APD90 << ", " << Cell.V_90 << ", " << Cell.APD_90 << ", " << Cell.DI << ", " << S1_cycle<< ", "<< Cell.peak_slope<< ", " << Cell.t_I_Na_peak << ", " << Cell.I_Na_peak << ", " << Cell.t_I_Na_L_peak << ", " << Cell.I_Na_L_peak <<  std::endl;
                
                        }
            
                        if (((t>=0.80*waitTime) && (t<0.80*waitTime+dt)) || ((t>=waitTime-1.0) && (t<waitTime-1.0+dt))) {//Last data    collection is at waitTime - 1 to ensure the stimulus current has not been added yet.
                            result_file3 << Cell.V << ", " << Cell.dV << ", " << Cell.Na_in<< ", " << Cell.K_in << ", " << Cell.Ca_in << ", " << Cell.Ca_sr << ", "<< Cell.Ca_ss << ", " << Cell.Ca_in_buffer << ", "<< Cell.Ca_ss_buffer << ", " << Cell.Ca_sr_buffer << ", "
                            << Cell.I_Na << ", "<<Cell.I_Na_L<< ", " << Cell.I_Ca_L << ", " << Cell.I_Kr << ", " << Cell.I_Ks << ", " << Cell.I_K1<< ", "<< Cell.I_Kp << ", " << Cell.I_to << ", " << Cell.I_Na_Ca << ", " << Cell.I_Na_K << ", " << Cell.I_p_Ca << ", " << Cell.I_Ca_b << ", " << Cell.I_Na_b << ", " << Cell.I_Na_ion_total << ", " << Cell.I_Ca_ion_total << ", " << Cell.I_K_ion_total << ", " << Cell.I_total << ", " << Cell.I_tr << ", " << Cell.I_leak << ", " << Cell.I_up << ", " << Cell.I_rel << ", " << Cell.m << ", " << Cell.h << ", " << Cell.j << ", " << Cell.d << ", " << Cell.f << ", " << Cell.f2 << ", " << Cell.f_Ca << ", " << Cell.r << ", " << Cell.s << ", " << Cell.xr1 << ", " << Cell.xr2 << ", " << Cell.xs << ", " << Cell.OO << ", " << Cell.R_bar << std::endl;
                
                        }
            
                        if ((counter%10==0) && (S1_cycle <=50) && (S1_cycle > 0)){result_file4 << t << ", " << Cell.V << ", " << Cell.I_Na<< ", " << Cell.I_Ca_L << ", " << Cell.I_Na_L << ", "<< Cell.I_up << ", "<< Cell.I_leak << ", "  << Cell.I_Ks << ", " << Cell.Ca_in << ", " << Cell.Ca_ss << std::endl;}// End of individual result files
            
            
                        /*******************************************************************************************************************************************************
                         End of result files
                         *******************************************************************************************************************************************************/
            
                        counter++;
            
                        if (fmod(t,5000)<dt) {std::cout << "S1_cycle=" <<S1_cycle<< ", " << "t=" << t << ", " "interval =" << interval  << std::endl;}
                
                    }//End of For tTemp ==
                    tTemp = 0;
        
                }//End of for S1
                S1_cycle = 1;    //Resets S1 cycle number
    
                //}    //End of for interval
        
                result_file.close();
                result_file0.close();
                result_file3.close();
                result_file4.close();
                
        }//end of loop for drug concentration
    }
        
    }
    
}

#endif /* Na_current_Dynamics_h */
