//
//  WT_Lido_Implicit.h
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in "Rate-dependent effects of lidocaine on cardiac dynamics: Development and analysis of a low-dimensional drug-channel interaction model" PLOS Computational Biology
//

////This code for the Ten Tusscher model was originally copied from code attached to "Jonathan D. Moreno, Z. Iris Zhu, Pei-Chi Yang, John R. Bankston, Mao-Tsuen Jeng, Chaoyi Kang, Lian- guo Wang, Jason D. Bayer, David J. Christini, Natalia A. Trayanova, Crystal M. Ripplinger, Robert S. Kass, and Colleen E. Clancy. A computational model to predict the effects of class 1 anti-arrhythmic drugs on ventricular rhythms. Science Translational Medicine, 3, 2011."

#ifndef WT_Lido_Implicit_h
#define WT_Lido_Implicit_h

#include <iostream>
#include <math.h>
#include "Global_variables.h"

// #include <Accelerate/Accelerate.h>
// #include "mat_inv.h"
// #include "solve_LU.h"

// extern "C" void dgetrf( int * M, int * N, double * A, int * LDA, int * IPIV, int * info);
// extern "C" void dgetri( int * N, double * A, int * LDA, int * IPIV, double * Work, int * LWork, int *info);

using namespace std;


void WT_SCN5A_Lidocaine_function(Cell_param *Cell_ptr, double t, double k_on_fac, double k_off_fac);

void WT_SCN5A_Lidocaine_function(Cell_param *Cell_ptr, double t, double k_on_fac, double k_off_fac){
    
    double E_Na, sum;
    
    double O , OS , C1 , C2 , C3 , IC3 , IC2 , IF , IM1 , IM2 ;
    double DO , DOS , DC1 , DC2 , DC3 , DIC3 , DIC2 , DIF , DIM1 , DIM2 ;
    double D_O , D_OS , D_C1 , D_C2 , D_C3 , D_IC3 , D_IC2 , D_IF , D_IM1 , D_IM2 ;
    
    double O_n , OS_n , C1_n , C2_n , C3_n , IC3_n , IC2_n , IF_n , IM1_n , IM2_n ;
    double DO_n , DOS_n , DC1_n , DC2_n , DC3_n , DIC3_n , DIC2_n , DIF_n , DIM1_n , DIM2_n ;
    double D_O_n , D_OS_n , D_C1_n , D_C2_n , D_C3_n , D_IC3_n , D_IC2_n , D_IF_n , D_IM1_n , D_IM2_n ;
    
    double O_o , OS_o , C1_o , C2_o , C3_o , IC3_o , IC2_o , IF_o , IM1_o , IM2_o ;
    double DO_o , DOS_o , DC1_o , DC2_o , DC3_o , DIC3_o , DIC2_o , DIF_o , DIM1_o , DIM2_o ;
    double D_O_o , D_OS_o , D_C1_o , D_C2_o , D_C3_o , D_IC3_o , D_IC2_o , D_IF_o , D_IM1_o , D_IM2_o ;
    
    //    double coef_O , coef_OS , coef_C1 , coef_C2 , coef_C3 , coef_IC3 , coef_IC2 , coef_IF , coef_IM1 , coef_IM2 ;
    //    double coef_DO , coef_DOS , coef_DC1 , coef_DC2 , coef_DC3 , coef_DIC3 , coef_DIC2 , coef_DIF , coef_DIM1 , coef_DIM2 ;
    //    double coef_D_O , coef_D_OS , coef_D_C1 , coef_D_C2 , coef_D_C3 , coef_D_IC3 , coef_D_IC2 , coef_D_IF , coef_D_IM1 , coef_D_IM2 ;
    
    //    double MB[30][30], MF[30][30], MT[30][30], invMB[30][30], nextVect[30], Vect[30];
    //    int id1, id2, id3;
    
    // For Mac using Accelerate.h
    //    __CLPK_doublereal x,y;
    //    __CLPK_integer info, LWork;
    //    __CLPK_integer M = 30;
    //    __CLPK_integer N = 30;
    //    __CLPK_integer LDA = 30;
    //    __CLPK_integer IPIV[30];
    //    __CLPK_doublereal A[900], Work[30];
    
    //    double x,y;
    //    int info, LWork;
    //    int M = 30;
    //    int N = 30;
    //    int LDA = 30;
    //    int IPIV[30];
    //    double A[900], Work[30], invA[900];
    
    
    double a11, a12, a13, a2, a3, a4, a5;
    double b11, b12, b13, b2, b3, b4, b5;
    double ax, bx, ax1, bx1, ax2, bx2;
    double a13c, b13c, a13n, b13n;
    double a22, b22, a33, b33, a44, b44, a55, b55;
    double a_22, b_22, a_33, b_33, a_44, b_44, a_55, b_55;
    
    double kd, kon, k_on, koff, k_off, kcon, kc_on, kcoff, kc_off, ki_on, ki_off;
    double kd_closed, kd_open;
    double hh;
    
    hh = dt;
    const double G_Na=15.0;    //27.5  //18.5
    
    const double Q10=3;
    
    const double Tfactor=1.0/(pow(Q10, (37.0-(T-273))/10.0));
    
    const double pH=7.4;
    const double pKa=7.6;
    const double portion = 1/(1+ pow(10, (pH-pKa)) );
    const double diffusion=500;
    
    //const double Drug=0*(1e-6);
    //const double Drug = Cell_ptr->Drug;
    
    const double Drug_charged=Drug*portion;
    const double Drug_neutral=Drug*(1-portion);
    const double dd= -0.7;
    
    E_Na = (R*T/F)*log(Na_out/Cell_ptr->Na_in);
    
    
    //Rate Constants *************************************************************************************************************************************************************************
    
    //WT Fits Reduced Model (no IM1, IM2)
    a11= Tfactor*8.5539/(7.4392e-2*exp(-Cell_ptr->V/17.0)+ 2.0373e-1*exp(-Cell_ptr->V/150));
    a12= Tfactor*8.5539/(7.4392e-2*exp(-Cell_ptr->V/15.0)+ 2.0373e-1*exp(-Cell_ptr->V/150));
    a13= Tfactor*8.5539/(7.4392e-2*exp(-Cell_ptr->V/12.0)+ 2.0373e-1*exp(-Cell_ptr->V/150));
    b11= Tfactor*7.5215e-2*exp(-Cell_ptr->V/20.3);
    b12= Tfactor*2.7574*exp(-(Cell_ptr->V-5)/20.3);
    b13= Tfactor*4.7755e-1*exp(-(Cell_ptr->V-10)/20.3);
    
    a3 = Tfactor*5.1458e-6*exp(-Cell_ptr->V/8.2471);
    b3=Tfactor*6.1205*exp((Cell_ptr->V)/13.542);
    
    
    a2= Tfactor*(13.370*exp(Cell_ptr->V/43.749));
    b2= ((a13*a2*a3)/(b13*b3));
    
    a4 = 0*a2;
    b4 = 0*a3;
    a5= 0*a2;
    b5 = 0*a3;
    
    ax = 3.4229e-2*a2;
    bx = 1.7898e-2*a3;
    
    
    ax1 = 6.3992e-07  *ax;
    bx1 = 1.3511e+00  *bx;
    a13c = 5.6974e-03  *a13;
    a22 = 6.7067e-06  *a2;
    b33 =  1.9698e-05 *b3;
    a33 = 3.2976e+00  *a3;
    
    a44 = 0  *a2;
    b44 = 0  *a3;
    a55 = b55 = 0;
    
    ax2 = 1.3110e-01  *ax;
    a13n = 8.4559e+01  *a13;
    a_22 = 1.7084e-05  *a2;
    b_33 = 4.8477e+00  *b3;
    
    a_44 = 0  *a2;
    b_44 = 0  *a3;
    a_55 = b_55 = 0;
    
    
    
    const double kd0=318*(1e-6);
    kd_open=kd0*exp( (dd*Cell_ptr->V*F) /(R*T));
    
    
    // charged drug
    kon=Drug_charged*diffusion*k_on_fac;
    koff=kd_open*diffusion*k_off_fac;
    kcoff = koff*k_off_fac;
    kcon = kon*k_on_fac;
    
    if (Drug ==0 || Drug_charged ==0 ){b13c = 0;}
    else{b13c = (b13*kcon*koff*a13c)/(kon*kcoff*a13);}
    
    if (b13c ==0){b22 = 0;}
    else {b22=(a13c*a22*a33)/(b13c*b33);}
    
    // neutral drug
    k_on = Drug_neutral*diffusion*k_on_fac;
    k_off=400*(1e-6)*diffusion*k_off_fac;
    ki_on=k_on/2*k_on_fac;
    ki_off=3.4*(1e-6)*diffusion*k_off_fac;
    kc_on=k_on/2*k_on_fac;
    kc_off=900*(1e-6)*diffusion*k_off_fac;
    
    if (Drug ==0 || Drug_neutral ==0 ){a_33 = 0;}
    else {a_33 = (ki_off*a3*kc_on*b_33)/(ki_on*kc_off*b3);}
    
    if (Drug ==0 || Drug_neutral ==0){b13n = 0;}
    else {b13n = (b13*kc_on*a13n*k_off)/(kc_off*a13*k_on);}
    
    if (b13n==0){b_22 =0;}
    else {b_22 = (a_33*a13n*a_22)/(b_33*b13n);}
    
    if (Drug ==0 || Drug_neutral ==0){bx2 = 0;}
    else{bx2 = (bx*k_on*ax2*ki_off)/(ax*ki_on*k_off);}
    
    //Initial Conditions *************************************************************************************************************************************************************************
    if (t==0){
        
        //Drug Free States
        Cell_ptr->IC3 = Cell_ptr->IC3;
        Cell_ptr->IC2 = Cell_ptr->IC2;
        Cell_ptr->IF = Cell_ptr->IF;
        Cell_ptr->IM1 = Cell_ptr->IM1;
        Cell_ptr->IM2 = Cell_ptr->IM2;
        Cell_ptr->C3 = Cell_ptr->C3;
        Cell_ptr->C2 = Cell_ptr->C2;
        Cell_ptr->C1 = Cell_ptr->C1;
        Cell_ptr->O = Cell_ptr->O;
        Cell_ptr->OS =Cell_ptr->OS;
        //Charged Drug States
        Cell_ptr->DIC3 = Cell_ptr->DIC3;
        Cell_ptr->DIC2 =Cell_ptr->DIC2;
        Cell_ptr->DIF = Cell_ptr->DIF;
        Cell_ptr->DIM1 = Cell_ptr->DIM1;
        Cell_ptr->DIM2 = Cell_ptr->DIM2;
        Cell_ptr->DC3 = Cell_ptr->DC3;
        Cell_ptr->DC2 =Cell_ptr->DC2;
        Cell_ptr->DO = Cell_ptr->DO;
        Cell_ptr->DOS = Cell_ptr->DOS;
        Cell_ptr->DC1 = Cell_ptr->DC1;
        //Neutral Drug States
        Cell_ptr->D_IC3 = Cell_ptr->D_IC3;
        Cell_ptr->D_IC2 = Cell_ptr->D_IC2;
        Cell_ptr->D_IF = Cell_ptr->D_IF;
        Cell_ptr->D_IM1 = Cell_ptr->D_IM1;
        Cell_ptr->D_IM2 = Cell_ptr->D_IM2;
        Cell_ptr->D_C3 = Cell_ptr->D_C3;
        Cell_ptr->D_C2 = Cell_ptr->D_C2;
        Cell_ptr->D_O = Cell_ptr->D_O;
        Cell_ptr->D_OS = Cell_ptr->D_OS;
        Cell_ptr->D_C1 = Cell_ptr->D_C1;
        /*
         //Drug Free States
         Cell_ptr->IC3 = 0;
         Cell_ptr->IC2 = 0;
         Cell_ptr->IF = 0;
         Cell_ptr->IM1 = 0;
         Cell_ptr->IM2 = 0;
         Cell_ptr->C3 = 1;
         Cell_ptr->C2 = 0;
         Cell_ptr->O = 0;
         Cell_ptr->OS =0;
         //Charged Drug States
         Cell_ptr->DIC3 = 0;
         Cell_ptr->DIC2 =0;
         Cell_ptr->DIF = 0;
         Cell_ptr->DIM1 = 0;
         Cell_ptr->DIM2 = 0;
         Cell_ptr->DC3 = 0;
         Cell_ptr->DC2 =0;
         Cell_ptr->DO = 0;
         Cell_ptr->DOS = 0;
         Cell_ptr->DC1 = 0;
         //Neutral Drug States
         Cell_ptr->D_IC3 = 0;
         Cell_ptr->D_IC2 = 0;
         Cell_ptr->D_IF = 0;
         Cell_ptr->D_IM1 = 0;
         Cell_ptr->D_IM2 = 0;
         Cell_ptr->D_C3 = 0;
         Cell_ptr->D_C2 = 0;
         Cell_ptr->D_O = 0;
         Cell_ptr->D_OS = 0;
         Cell_ptr->D_C1 = 0;
         
         
         Cell_ptr->C1 = 1 - ( Cell_ptr->O + Cell_ptr->OS + Cell_ptr->C3 + Cell_ptr->C2 + Cell_ptr->IC3
         + Cell_ptr->IC2 + Cell_ptr->IF + Cell_ptr->IM1 + Cell_ptr->IM2 +
         Cell_ptr->DO + Cell_ptr->DOS + Cell_ptr->DC1 + Cell_ptr->DC2 + Cell_ptr->DC3
         + Cell_ptr->DIC3 + Cell_ptr->DIC2 + Cell_ptr->DIF + Cell_ptr->DIM1 + Cell_ptr->DIM2 +
         Cell_ptr->D_O + Cell_ptr->D_OS + Cell_ptr->D_C1 + Cell_ptr->D_C2 + Cell_ptr->D_C3
         + Cell_ptr->D_IC3 + Cell_ptr->D_IC2 + Cell_ptr->D_IF + Cell_ptr->D_IM1 + Cell_ptr->D_IM2);
         */
    } //End if t == O
    
    else if (t>0.0){
        
        IC3 = Cell_ptr->IC3;
        IC2 = Cell_ptr->IC2;
        IF = Cell_ptr->IF;
        IM1 = Cell_ptr->IM1;
        IM2 = Cell_ptr->IM2;
        C3 = Cell_ptr->C3;
        C2 = Cell_ptr->C2;
        C1 = Cell_ptr->C1;
        O = Cell_ptr->O;
        OS = Cell_ptr->OS;
        DC3 = Cell_ptr->DC3;
        DC2 = Cell_ptr->DC2;
        DC1 = Cell_ptr->DC1;
        DO = Cell_ptr->DO;
        DOS = Cell_ptr->DOS;
        DIC3 = Cell_ptr->DIC3;
        DIC2 = Cell_ptr->DIC2;
        DIF = Cell_ptr->DIF;
        DIM1 = Cell_ptr->DIM1;
        DIM2 = Cell_ptr->DIM2;
        D_C3 = Cell_ptr->D_C3;
        D_C2 = Cell_ptr->D_C2;
        D_C1 = Cell_ptr->D_C1;
        D_O = Cell_ptr->D_O;
        D_OS = Cell_ptr->D_OS;
        D_IC3 = Cell_ptr->D_IC3;
        D_IC2 = Cell_ptr->D_IC2;
        D_IF = Cell_ptr->D_IF;
        D_IM1 = Cell_ptr->D_IM1;
        D_IM2 = Cell_ptr->D_IM2;
        
        // Backward method
        // MB * nextVect_(n+1) = Vect_n
        // where MB is a matrix, and Vect_i is a vector
        
        // Set up the matrix
        //  MB = zeros(30);
        //        for (id1 = 0; id1 < 30; id1 += 1){
        //            for(id2 = 0; id2 < 30; id2 += 1) {
        //                MB[id1][id2] = 0;
        //            } // end of id2 looping
        //        } // end of id1 looping
        
        
        //        // Calculation of k parameters for the Runge Kutta interations
        //
        //        //Calculation of K1 *************************************************************************************************************************************************************************
        //        //Drug Free States
        //        O = hh*( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS  - O * (b13 + a2  + kon + k_on + ax));
        //        C1 = hh*(a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 - C1*(b12 + b3 + a13  + kcon + kc_on));
        //        C2 = hh*(a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 - C2*(b11 + b3 + a12  + kcon + kc_on));
        //        C3 = hh*(a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3*(b3 + a11  + kcon + kc_on));
        //        IC3 = hh*(b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3*(a11 + a3 + ki_on));
        //        IC2 = hh*(a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2*(b11 + a3 + a12 + ki_on));
        //        IF = hh*(a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF*(b12 + b2 + a3 + a4 + ki_on));
        //        IM1 = hh*(a4 * IF + b5 * IM2 - IM1*(b4 + a5));
        //        IM2 = hh*(a5 * IM1 - IM2*(b5));
        //        OS = hh*(ax * O + ki_off * D_OS - OS*(bx + ki_on));
        //
        //        //Charged Drug Bound States
        //        DO = hh*(kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO*(koff + b13c + a22 + ax1 ));
        //        DC1 = hh*(kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1*(kcoff + b12 + b33 + a13c ));
        //        DC2 = hh*(kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2*(kcoff + b11 + b33 + a12 ));
        //        DC3 = hh*(kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3*(kcoff+ b33 + a11 ));
        //        DOS = hh*(ax1 * DO - bx1 * DOS);
        //        DIC3 = hh*(b33 * DC3 + b11 * DIC2 - DIC3*(a11 + a33));
        //        DIC2 = hh*(b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2*(a33 + b11 + a12));
        //        DIF = hh*(b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF*(a33 + b12 + a44 + b22));
        //        DIM1 = hh*(a44 * DIF + b55 * DIM2 - DIM1*( b44 + a55));
        //        DIM2 = hh*(a55 * DIM1 - b55 * DIM2);
        //
        //        //Neutral Drug Bound States
        //        D_O = hh*(k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O*(k_off + b13n + a_22 + ax2 ));
        //        D_C1 = hh*(kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1*(kc_off + b12 + b_33 + a13n ));
        //        D_C2 = hh*(kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2*(kc_off + b11 + b_33 + a12 ));
        //        D_C3 = hh*(kc_on * C3 + a_33 * D_IC3 + b11 * D_C2  - D_C3*(kc_off + b_33 + a11 ));
        //        D_OS = hh*(ax2 * D_O + ki_on * OS - D_OS*(bx2 + ki_off));
        //        D_IC3 = hh*(b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3*(a_33 + a11 + ki_off));
        //        D_IC2 = hh*(b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2*(a_33 + b11 + a12 + ki_off));
        //        D_IF = hh*(b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF*(a_33 + a_44 + b_22 + b12 + ki_off));
        //        D_IM1 = hh*(a_44 * D_IF + b_55 * D_IM2 - D_IM1*(b_44 + a_55));
        //        D_IM2 = hh*(a_55 * D_IM1 - b_55 * D_IM2);
        
        const double coef_O = (b13 + a2  + kon + k_on + ax);
        const double coef_C1 = (b12 + b3 + a13  + kcon + kc_on);
        const double coef_C2 = (b11 + b3 + a12  + kcon + kc_on);
        const double coef_C3 = (b3 + a11  + kcon + kc_on);
        const double coef_IC3 = (a11 + a3 + ki_on);
        const double coef_IC2 = (b11 + a3 + a12 + ki_on);
        const double coef_IF = (b12 + b2 + a3 + a4 + ki_on);
        const double coef_IM1 = (b4 + a5);
        const double coef_IM2 = b5;
        const double coef_OS = (bx + ki_on);
        
        const double coef_DO = (koff + b13c + a22 + ax1 );
        const double coef_DC1 = (kcoff + b12 + b33 + a13c );
        const double coef_DC2 = (kcoff + b11 + b33 + a12 );
        const double coef_DC3 = (kcoff+ b33 + a11 );
        const double coef_DOS = bx1;
        const double coef_DIC3 = (a11 + a33);
        const double coef_DIC2 = (a33 + b11 + a12);
        const double coef_DIF = (a33 + b12 + a44 + b22);
        const double coef_DIM1 = ( b44 + a55 );
        const double coef_DIM2 = b55 ;
        
        const double coef_D_O = (k_off + b13n + a_22 + ax2 );
        const double coef_D_C1 = (kc_off + b12 + b_33 + a13n );
        const double coef_D_C2 = (kc_off + b11 + b_33 + a12 );
        const double coef_D_C3 = (kc_off + b_33 + a11 );
        const double coef_D_OS = (bx2 + ki_off);
        const double coef_D_IC3 = (a_33 + a11 + ki_off);
        const double coef_D_IC2 = (a_33 + b11 + a12 + ki_off);
        const double coef_D_IF = (a_33 + a_44 + b_22 + b12 + ki_off);
        const double coef_D_IM1 = (b_44 + a_55);
        const double coef_D_IM2 = b_55;
        
        double hh2 = hh/2;
        
        const double co_O = 1. / ( 1 + hh2 * ( b13 + a2  + kon + k_on + ax ) );
        const double co_C1 = 1. / ( 1 + hh2 * ( b12 + b3 + a13  + kcon + kc_on));
        const double co_C2 = 1. / ( 1 + hh2 * ( b11 + b3 + a12  + kcon + kc_on));
        const double co_C3 = 1. / ( 1 + hh2 * ( b3 + a11  + kcon + kc_on));
        const double co_IC3 = 1. / ( 1 + hh2 * ( a11 + a3 + ki_on));
        const double co_IC2 = 1. / ( 1 + hh2 * ( b11 + a3 + a12 + ki_on));
        const double co_IF = 1. / ( 1 + hh2 * ( b12 + b2 + a3 + a4 + ki_on));
        const double co_IM1 = 1. / ( 1 + hh2 * ( b4 + a5));
        const double co_IM2 = 1. / ( 1 + hh2 * ( b5 ));
        const double co_OS = 1. / ( 1 + hh2 * ( bx + ki_on));
        
        const double co_DO = 1. / ( 1 + hh2 * ( koff + b13c + a22 + ax1 ));
        const double co_DC1 = 1. / ( 1 + hh2 * ( kcoff + b12 + b33 + a13c ));
        const double co_DC2 = 1. / ( 1 + hh2 * ( kcoff + b11 + b33 + a12 ));
        const double co_DC3 = 1. / ( 1 + hh2 * ( kcoff+ b33 + a11 ));
        const double co_DOS = 1. / ( 1 + hh2 * ( bx1));
        const double co_DIC3 = 1. / ( 1 + hh2 * ( a11 + a33));
        const double co_DIC2 = 1. / ( 1 + hh2 * ( a33 + b11 + a12));
        const double co_DIF = 1. / ( 1 + hh2 * ( a33 + b12 + a44 + b22));
        const double co_DIM1 = 1. / ( 1 + hh2 * (  b44 + a55 ));
        const double co_DIM2 = 1. / ( 1 + hh2 * ( b55 ) );
        
        const double co_D_O = 1. / ( 1 + hh2 * ( k_off + b13n + a_22 + ax2 ));
        const double co_D_C1 = 1. / ( 1 + hh2 * ( kc_off + b12 + b_33 + a13n ));
        const double co_D_C2 = 1. / ( 1 + hh2 * ( kc_off + b11 + b_33 + a12 ));
        const double co_D_C3 = 1. / ( 1 + hh2 * ( kc_off + b_33 + a11 ));
        const double co_D_OS = 1. / ( 1 + hh2 * ( bx2 + ki_off));
        const double co_D_IC3 = 1. / ( 1 + hh2 * ( a_33 + a11 + ki_off));
        const double co_D_IC2 = 1. / ( 1 + hh2 * ( a_33 + b11 + a12 + ki_off));
        const double co_D_IF = 1. / ( 1 + hh2 * ( a_33 + a_44 + b_22 + b12 + ki_off));
        const double co_D_IM1 = 1. / ( 1 + hh2 * ( b_44 + a_55));
        const double co_D_IM2 = 1. / ( 1 + hh2 * ( b_55 ));
        
        
        //Drug Free States
        O_o = O + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS  - O * coef_O );
        C1_o = C1 + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 - C1 * coef_C1 );
        C2_o = C2 + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 - C2 * coef_C2 );
        C3_o = C3 + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 - C3 * coef_C3 );
        IC3_o = IC3 + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3 - IC3 * coef_IC3 );
        IC2_o = IC2 + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2 - IC2 * coef_IC2 );
        IF_o = IF + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF - IF * coef_IF );
        IM1_o = IM1 + hh2 * (a4 * IF + b5 * IM2 - IM1 * coef_IM1 );
        IM2_o = IM2 + hh2 * (a5 * IM1 - IM2 * coef_IM2 );
        OS_o = OS + hh2 * (ax * O + ki_off * D_OS - OS * coef_OS );
        
        //Charged Drug Bound States
        DO_o = DO + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF  - DO * coef_DO );
        DC1_o = DC1 + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  - DC1 * coef_DC1 );
        DC2_o = DC2 + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1  - DC2 * coef_DC2 );
        DC3_o = DC3 + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3  - DC3 * coef_DC3 );
        DOS_o = DOS + hh2 * (ax1 * DO -  DOS * coef_DOS );
        DIC3_o = DIC3 + hh2 * (b33 * DC3 + b11 * DIC2 - DIC3 * coef_DIC3 );
        DIC2_o = DIC2 + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF - DIC2 * coef_DIC2 );
        DIF_o = DIF + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO - DIF * coef_DIF );
        DIM1_o = DIM1 + hh2 * (a44 * DIF + b55 * DIM2 - DIM1 * coef_DIM1 );
        DIM2_o = DIM2 + hh2 * (a55 * DIM1 - DIM2 * coef_DIM2 );
        
        //Neutral Drug Bound States
        D_O_o = D_O + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS  - D_O * coef_D_O );
        D_C1_o = D_C1 + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O  - D_C1 * coef_D_C1 );
        D_C2_o = D_C2 + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  - D_C2 * coef_D_C2 );
        D_C3_o = D_C3 + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2  - D_C3 * coef_D_C3 );
        D_OS_o = D_OS + hh2 * (ax2 * D_O + ki_on * OS - D_OS * coef_D_OS );
        D_IC3_o = D_IC3 + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 - D_IC3 * coef_D_IC3 );
        D_IC2_o = D_IC2 + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2 - D_IC2 * coef_D_IC2 );
        D_IF_o = D_IF + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF - D_IF * coef_D_IF );
        D_IM1_o = D_IM1 + hh2 * (a_44 * D_IF + b_55 * D_IM2 - D_IM1 * coef_D_IM1 );
        D_IM2_o = D_IM2 + hh2 * (a_55 * D_IM1 - D_IM2 * coef_D_IM2 );
        
        int iter = 0;
        double err_sum = 1;
        while ( err_sum > 1E-100 && iter < 100 ) {
            
            //            //Drug Free States
            //            O_n = co_O * ( O_o - hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS ) );
            //            C1_n = co_C1 * ( C1_o - hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 ) );
            //            C2_n = co_C2 * ( C2_o - hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 ) );
            //            C3_n = co_C3 * ( C3_o - hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 ) );
            //            IC3_n = co_IC3 * ( IC3_o - hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3  ) );
            //            IC2_n = co_IC2 * ( IC2_o - hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2  ) );
            //            IF_n = co_IF * ( IF_o - hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF  ) );
            //            IM1_n = co_IM1 * ( IM1_o - hh2 * (a4 * IF + b5 * IM2  ) );
            //            IM2_n = co_IM2 * ( IM2_o - hh2 * (a5 * IM1  ) );
            //            OS_n = co_OS * ( OS_o - hh2 * (ax * O + ki_off * D_OS  ) );
            //
            //            //Charged Drug Bound States
            //            DO_n = co_DO * ( DO_o - hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF   ) );
            //            DC1_n = co_DC1 * ( DC1_o - hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  ) );
            //            DC2_n = co_DC2 * ( DC2_o - hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1   ) );
            //            DC3_n = co_DC3 * ( DC3_o - hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3   ) );
            //            DOS_n = co_DOS * ( DOS_o - hh2 * (ax1 * DO  ) );
            //            DIC3_n = co_DIC3 * ( DIC3_o - hh2 * (b33 * DC3 + b11 * DIC2 ) );
            //            DIC2_n = co_DIC2 * ( DIC2_o - hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF  ) );
            //            DIF_n = co_DIF * ( DIF_o - hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ) );
            //            DIM1_n = co_DIM1 * ( DIM1_o - hh2 * (a44 * DIF + b55 * DIM2  ) );
            //            DIM2_n = co_DIM2 * ( DIM2_o - hh2 * (a55 * DIM1 ) );
            //
            //            //Neutral Drug Bound States
            //            D_O_n = co_D_O * ( D_O_o - hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS   ) );
            //            D_C1_n = co_D_C1 * ( D_C1_o - hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O   ) );
            //            D_C2_n = co_D_C2 * ( D_C2_o - hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  ) );
            //            D_C3_n = co_D_C3 * ( D_C3_o - hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2   ) );
            //            D_OS_n = co_D_OS * ( D_OS_o - hh2 * (ax2 * D_O + ki_on * OS ) );
            //            D_IC3_n = co_D_IC3 * ( D_IC3_o - hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 ) );
            //            D_IC2_n = co_D_IC2 * ( D_IC2_o - hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2  ) );
            //            D_IF_n = co_D_IF * ( D_IF_o - hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF  ) );
            //            D_IM1_n = co_D_IM1 * ( D_IM1_o - hh2 * (a_44 * D_IF + b_55 * D_IM2  ) );
            //            D_IM2_n = co_D_IM2 * ( D_IM2_o - hh2 * (a_55 * D_IM1  ) );
            
            //Drug Free States
            O_n = co_O * ( O_o + hh2 * ( a13 * C1 + b2 * IF + koff * DO + k_off * D_O + bx * OS ) );
            C1_n = co_C1 * ( C1_o + hh2 * (a12 * C2 + a3 * IF + b13 * O  + kcoff * DC1 + kc_off * D_C1 ) );
            C2_n = co_C2 * ( C2_o + hh2 * (a11 * C3 + a3 * IC2 + b12 * C1  + kcoff * DC2 + kc_off * D_C2 ) );
            C3_n = co_C3 * ( C3_o + hh2 * (a3 * IC3 + b11 * C2  + kcoff * DC3 + kc_off * D_C3 ) );
            IC3_n = co_IC3 * ( IC3_o + hh2 * (b3 * C3 + b11 * IC2 + ki_off * D_IC3  ) );
            IC2_n = co_IC2 * ( IC2_o + hh2 * (a11 * IC3 + b3 * C2 + b12 * IF + ki_off * D_IC2  ) );
            IF_n = co_IF * ( IF_o + hh2 * (a12 * IC2 + b3 * C1 + b4 * IM1 + a2 * O + ki_off * D_IF  ) );
            IM1_n = co_IM1 * ( IM1_o + hh2 * (a4 * IF + b5 * IM2  ) );
            IM2_n = co_IM2 * ( IM2_o + hh2 * (a5 * IM1  ) );
            OS_n = co_OS * ( OS_o + hh2 * (ax * O + ki_off * D_OS  ) );
            
            //Charged Drug Bound States
            DO_n = co_DO * ( DO_o + hh2 * (kon * O + a13c * DC1 + bx1 * DOS + b22 * DIF   ) );
            DC1_n = co_DC1 * ( DC1_o + hh2 * (kcon * C1 + a12 * DC2 + a33 * DIF + b13c * DO  ) );
            DC2_n = co_DC2 * ( DC2_o + hh2 * (kcon * C2 + a11 * DC3 + a33 * DIC2 + b12 * DC1   ) );
            DC3_n = co_DC3 * ( DC3_o + hh2 * (kcon * C3 + b11 * DC2 + a33 * DIC3   ) );
            DOS_n = co_DOS * ( DOS_o + hh2 * (ax1 * DO  ) );
            DIC3_n = co_DIC3 * ( DIC3_o + hh2 * (b33 * DC3 + b11 * DIC2 ) );
            DIC2_n = co_DIC2 * ( DIC2_o + hh2 * (b33 * DC2 + a11 * DIC3 + b12 * DIF  ) );
            DIF_n = co_DIF * ( DIF_o + hh2 * (b33 * DC1 + a12 * DIC2 + b44 * DIM1 + a22 * DO ) );
            DIM1_n = co_DIM1 * ( DIM1_o + hh2 * (a44 * DIF + b55 * DIM2  ) );
            DIM2_n = co_DIM2 * ( DIM2_o + hh2 * (a55 * DIM1 ) );
            
            //Neutral Drug Bound States
            D_O_n = co_D_O * ( D_O_o + hh2 * (k_on * O + a13n * D_C1 + b_22 * D_IF + bx2 * D_OS   ) );
            D_C1_n = co_D_C1 * ( D_C1_o + hh2 * (kc_on * C1 + a12 * D_C2 + a_33 * D_IF + b13n * D_O   ) );
            D_C2_n = co_D_C2 * ( D_C2_o + hh2 * (kc_on * C2 + a11 * D_C3 + a_33 * D_IC2 + b12 * D_C1  ) );
            D_C3_n = co_D_C3 * ( D_C3_o + hh2 * (kc_on * C3 + a_33 * D_IC3 + b11 * D_C2   ) );
            D_OS_n = co_D_OS * ( D_OS_o + hh2 * (ax2 * D_O + ki_on * OS ) );
            D_IC3_n = co_D_IC3 * ( D_IC3_o + hh2 * (b_33 * D_C3 + b11 * D_IC2 + ki_on * IC3 ) );
            D_IC2_n = co_D_IC2 * ( D_IC2_o + hh2 * (b_33 * D_C2 + a11 * D_IC3 + b12 * D_IF + ki_on * IC2  ) );
            D_IF_n = co_D_IF * ( D_IF_o + hh2 * (b_33 * D_C1 + a12 * D_IC2 + b_44 * D_IM1 + a_22 * D_O + ki_on * IF  ) );
            D_IM1_n = co_D_IM1 * ( D_IM1_o + hh2 * (a_44 * D_IF + b_55 * D_IM2  ) );
            D_IM2_n = co_D_IM2 * ( D_IM2_o + hh2 * (a_55 * D_IM1  ) );
            
            
            err_sum = fabs( IC3 - IC3_n ) + fabs( IC2 - IC2_n ) +  + fabs( IF - IF_n ) +  + fabs( IM1 - IM1_n ) +  + fabs( IM2 - IM2_n ) +  + fabs( C3 - C3_n ) +  + fabs( C2 - C2_n ) +  + fabs( C1 - C1_n ) +  + fabs( O - O_n ) +  + fabs( OS - OS_n ) +  + fabs( DC3 - DC3_n ) +  + fabs( DC2 - DC2_n ) +  + fabs( DC1 - DC1_n ) +  + fabs( DO - DO_n ) +  + fabs( DOS - DOS_n ) +   + fabs( DIC3 - DIC3_n ) +  + fabs( DIC2 - DIC2_n ) +  + fabs( DIF - DIF_n ) +  + fabs( DIM1 - DIM1_n ) +  + fabs( DIM2 - DIM2_n ) +  + fabs( D_C3 - D_C3_n ) +  + fabs( D_C2 - D_C2_n ) +  + fabs( D_C1 - D_C1_n ) +  + fabs( D_O - D_O_n ) +  + fabs( D_OS - D_OS_n ) +  + fabs( D_IC3 - D_IC3_n ) +  + fabs( D_IC2 - D_IC2_n ) +  + fabs( D_IF - D_IF_n ) +  + fabs( D_IM1 - D_IM1_n ) +  + fabs( D_IM2 - D_IM2_n );
            
            
            IC3 = IC3_n;
            IC2 = IC2_n;
            IF = IF_n;
            IM1 = IM1_n;
            IM2 = IM2_n;
            C3 = C3_n;
            C2 = C2_n;
            C1 = C1_n;
            O = O_n;
            OS = OS_n;
            DC3 = DC3_n;
            DC2 = DC2_n;
            DC1 = DC1_n;
            DO = DO_n;
            DOS = DOS_n;
            DIC3 = DIC3_n;
            DIC2 = DIC2_n;
            DIF = DIF_n;
            DIM1 = DIM1_n;
            DIM2 = DIM2_n;
            D_C3 = D_C3_n;
            D_C2 = D_C2_n;
            D_C1 = D_C1_n;
            D_O = D_O_n;
            D_OS = D_OS_n;
            D_IC3 = D_IC3_n;
            D_IC2 = D_IC2_n;
            D_IF = D_IF_n;
            D_IM1 = D_IM1_n;
            D_IM2 = D_IM2_n;
            
            
            iter++;
        }
        
        // cout << iter << "\t" << err_sum << endl << endl;
        
        /******************************************************************************************************************************************************************************************************************/
        //        //        MF = eye(30) - 0.5 * dt * MB;
        //        //        MB = 0.5 * dt * MB + eye(30);
        //        for (id1 = 0; id1 < 30; id1 += 1){
        //            for (id2 = 0; id2 < 30; id2 += 1) {
        //                MB[id1][id2] = 0.5 * dt * MB[id1][id2];
        //                MF[id1][id2] = -MB[id1][id2];
        //            }
        //            MB[id1][id1] = MB[id1][id1] + 1;
        //            MF[id1][id1] = MF[id1][id1] + 1;
        //        }
        
        
        
        //        // Compute inv(MB) : Inverse of the matrix MB
        //        // Get LU factorization first,
        //        // then use LU to compute inv(MB)
        //        id3 = 0;
        //        for (id1 = 0; id1 < 30; id1 += 1) {
        //            for (id2 = 0; id2 < 30; id2 += 1){
        //                // For Mac using Accelerate.h
        //                // Fortran is Column major, C++ is Row major
        //                //  A[id3] = MB[id2][id1];
        //
        //                // For MTJ C++ subroutine
        //                A[id3] = MB[id1][id2];
        //
        //                id3 += 1;
        //            }
        //        }
        //
        //        LWork = N;
        //        // For Mac using Accelerate.h
        //        // dgetrf_( &M, &N, A, &LDA, IPIV, &info);
        //        // dgetri_( &N, A, &LDA, IPIV, Work, &LWork, &info);
        //
        //
        //        // dgetrf( &M, &N, A, &LDA, IPIV, &info);
        //        // dgetri( &N, A, &LDA, IPIV, Work, &LWork, &info);
        //
        //        // For MTJ C++ subroutine
        //        mat_inv( M, A, invA );
        //
        //
        //        // if (info != 0) { return info;}
        //        // Let invMB = inv(MB)
        //        id3 = 0;
        //        for (id1 = 0; id1 < 30; id1 += 1) {
        //            for (id2 = 0; id2 < 30; id2 += 1){
        //                // Fortran is Column major, C++ is Row major
        //                //  invMB[id2][id1] = A[id3];
        //
        //                // For MTJ C++ subroutine
        //                 invMB[id1][id2] = invA[id3];
        //
        //                id3 += 1;
        //                // cout << invMB[id1][id2] << "\t";
        //            }
        //            // cout << endl;
        //        }
        //
        //
        //
        //        for (id1 = 0; id1 < 30; id1 += 1) {
        //            for (id2 = 0; id2 < 30; id2 += 1) {
        //                MT[id1][id2] = 0;
        //                for(id3 = 0; id3 < 30; id3 += 1) {
        //                    MT[id1][id2] += invMB[id1][id3] * MF[id3][id2];
        //                }
        //            }
        //        }
        //
        //
        //
        //
        //        // Compute Vect_(n+1) = inv(MB) * Vect_n
        //        for (id1 = 0; id1 < 30; id1 += 1){
        //            nextVect[id1] = 0;
        //            for(id2 = 0; id2 < 30; id2 += 1) {
        //                // the value for the entry id2 of the vector
        //                nextVect[id1] += MT[id1][id2] * Vect[id2];
        //            } // end of id2 looping
        //        } // end of id1 looping
        //
        //        //#pragma omp critical
        //        //{
        //        //        for ( id1 = 0; id1 < 30; id1++ ) {
        //        //        cout << Vect[id1] << ", ";
        //        //        }
        //        //    cout << endl << endl;
        //        //}
        //
        //
        //        // Vect = [IC3; IC2; IF; IM1; IM2; C3; C2; C1; O; OS; DC3; DC2; DC1; DO; DOS; ...
        //        //     DIC3; DIC2; DIF; DIM1; DIM2; D_C3; D_C2; D_C1; D_O; D_OS; ...
        //        //     D_IC3; D_IC2; D_IF; D_IM1; D_IM2];
        //
        //        //        nextVect = inv(MB) * MF * Vect;
        
        
        //        Cell_ptr->IC3 = nextVect[0];
        //        Cell_ptr->IC2 = nextVect[1];
        //        Cell_ptr->IF = nextVect[2];
        //        Cell_ptr->IM1 = nextVect[3];
        //        Cell_ptr->IM2 = nextVect[4];
        //        Cell_ptr->C3 = nextVect[5];
        //        Cell_ptr->C2 = nextVect[6];
        //        Cell_ptr->C1 = nextVect[7];
        //        Cell_ptr->O = nextVect[8];
        //        Cell_ptr->OS = nextVect[9];
        //        Cell_ptr->DC3 = nextVect[10];
        //        Cell_ptr->DC2 = nextVect[11];
        //        Cell_ptr->DC1 = nextVect[12];
        //        Cell_ptr->DO = nextVect[13];
        //        Cell_ptr->DOS = nextVect[14];
        //        Cell_ptr->DIC3 = nextVect[15];
        //        Cell_ptr->DIC2 = nextVect[16];
        //        Cell_ptr->DIF = nextVect[17];
        //        Cell_ptr->DIM1 = nextVect[18];
        //        Cell_ptr->DIM2 = nextVect[19];
        //        Cell_ptr->D_C3 = nextVect[20];
        //        Cell_ptr->D_C2 = nextVect[21];
        //        Cell_ptr->D_C1 = nextVect[22];
        //        Cell_ptr->D_O = nextVect[23];
        //        Cell_ptr->D_OS = nextVect[24];
        //        Cell_ptr->D_IC3 = nextVect[25];
        //        Cell_ptr->D_IC2 = nextVect[26];
        //        Cell_ptr->D_IF = nextVect[27];
        //        Cell_ptr->D_IM1 = nextVect[28];
        //        Cell_ptr->D_IM2 = nextVect[29];
        
        Cell_ptr->IC3 = IC3_n;
        Cell_ptr->IC2 = IC2_n;
        Cell_ptr->IF = IF_n;
        Cell_ptr->IM1 = IM1_n;
        Cell_ptr->IM2 = IM2_n;
        Cell_ptr->C3 = C3_n;
        Cell_ptr->C2 = C2_n;
        Cell_ptr->C1 = C1_n;
        Cell_ptr->O = O_n;
        Cell_ptr->OS = OS_n;
        Cell_ptr->DC3 = DC3_n;
        Cell_ptr->DC2 = DC2_n;
        Cell_ptr->DC1 = DC1_n;
        Cell_ptr->DO = DO_n;
        Cell_ptr->DOS = DOS_n;
        Cell_ptr->DIC3 = DIC3_n;
        Cell_ptr->DIC2 = DIC2_n;
        Cell_ptr->DIF = DIF_n;
        Cell_ptr->DIM1 = DIM1_n;
        Cell_ptr->DIM2 = DIM2_n;
        Cell_ptr->D_C3 = D_C3_n;
        Cell_ptr->D_C2 = D_C2_n;
        Cell_ptr->D_C1 = D_C1_n;
        Cell_ptr->D_O = D_O_n;
        Cell_ptr->D_OS = D_OS_n;
        Cell_ptr->D_IC3 = D_IC3_n;
        Cell_ptr->D_IC2 = D_IC2_n;
        Cell_ptr->D_IF = D_IF_n;
        Cell_ptr->D_IM1 = D_IM1_n;
        Cell_ptr->D_IM2 = D_IM2_n;
        
        
        /******************************************************************************************************************************************************************************************************************/
        
        
        sum = Cell_ptr->O + Cell_ptr->OS + Cell_ptr->C1 + Cell_ptr->C2 + Cell_ptr->C3 + Cell_ptr->IC3 + Cell_ptr->IC2 + Cell_ptr->IF + Cell_ptr->IM1 + Cell_ptr->IM2 +
        Cell_ptr->DO + Cell_ptr->DOS + Cell_ptr->DC1 + Cell_ptr->DC2 + Cell_ptr->DC3 + Cell_ptr->DIC3 + Cell_ptr->DIC2 + Cell_ptr->DIF + Cell_ptr->DIM1 + Cell_ptr->DIM2 +
        Cell_ptr->D_O + Cell_ptr->D_OS + Cell_ptr->D_C1 + Cell_ptr->D_C2 + Cell_ptr->D_C3 + Cell_ptr->D_IC3 + Cell_ptr->D_IC2 + Cell_ptr->D_IF + Cell_ptr->D_IM1  + Cell_ptr->D_IM2;
        
        // cout << "sum = " << Cell_ptr->sum << endl;
        if (fabs (sum -1.0) > 0.0001) {cout << "Error in Wild Type Sum at t = " << t << ",\t sum = " << sum << endl;}
        
        
    } //End else if t>0
    
    Cell_ptr->I_Na = G_Na*(Cell_ptr->O)*(Cell_ptr->V - E_Na);
}

#endif /* WT_Lido_Implicit_h */
