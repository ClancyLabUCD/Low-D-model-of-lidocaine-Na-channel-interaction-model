//
//  main.cpp
//  ten_Tusscher_single_cell_with_Moreno_2011_and_Simple_Lidocaine_model
//
//  Created by Steffen Docken on 4/6/18.
//  Copyright © 2018 Steffen Docken (Lewis Lab). All rights reserved.
//  Published in "Rate-dependent effects of lidocaine on cardiac dynamics: Development and analysis of a low-dimensional drug-channel interaction model" PLOS Computational Biology
//

#include <iostream>
#include "Na_current_Dynamics.h"
#include "Upstroke_vs_period.h"
#include "Moreno_Drug_Dynamics.h"

int main() {
    //Na_current_Dynamics();
    Upstroke_vs_period();
    //Moreno_Drug_Dynamics();
    
    return 0;
}
