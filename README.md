# Low-D-model-of-lidocaine-Na-channel-interaction-model

C++ code for low-dimensional model of lidocaine-Na channel interaction model


The most important file here is I_Na_Simple_Lidocaine.h as this file contains the code for our low-dimensional lidocaine-Na channel interaction model.


Files for running simulations:

main.cpp: This file declares the protocol to be simulated by calling the appropriate function

Upstroke_vs_period.h: This code paces the ten Tusscher et al. 2004 and 2006 model with various Na current models (ten Tusscher et al. 2004 and 2006, our model, or Moreno et al. 2011) 500 times at various BCL to calculate the upstroke velocity as a function of BCL during steady pacing.

Na_current_Dynamics.h: This program brings the cell to steady state and then paces it 500 times to record dynamics of the Na-current variables and certain other cell variables over time.

Moreno_Drug_Dynamics.h: This program simulates a paced cell and records the fraction of Na channels in various conformational states and bound to neutral or charged lidocaine as predicted by the Moreno et al. 2011 model.


Files containing model equations and parameters:

Global_variables.h: Contains the global parameters of the ten Tusscher et al. 2004 and 2006 Human Ventricular Electrophysiology model along with our Na channel model and the Moreno et al. 2011 Na channel model.

Non_Na_currents.h: Contains the models for all ionic currents of ten Tusscher et al. 2004 and 2006 except the Na current model.

Numerical_Methods.h: Contains functions for updating the value of fast (relative to the time step) gating variables.

Ion_conc_dynamics.h: Contains the equations for updating ionic concentrations in the ten Tusscher et al. 2004 and 2006 model.

Na_currents.h: This code calls the desired Na current model (ten Tusscher et al. 2004 and 2006, our model, or Moreno et al. 2011) based on the value of the Na_model variable set in Global_variables.h.

WT_Lido_Implicit.h: This code runs the Moreno et al. lidocaine-Na channel interaction model.

Full_cell_functions.h: This file contains the functions that update Voltage and other full cell variables at the end of each iteration.

I_Na_TT.h: This code runs the ten Tusscher et al. 2004 and 2006 Na channel model.
