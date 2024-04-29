/*
    * @author Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 1-Dimensional Central Upwind Scheme
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Given the Conserved variables we can convert them to get the primitive variables

    ! Conserved Variables - Density(rho), Momentum(m), Energy(E)
    ! Primitive Variables - Density(rho), Velocity(u), Pressure(p)
*/
#pragma once

#include<iostream>
#include<vector>

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace PrimitiveVariables{
    vector< vector<double> > get_primitive_variables(vector< vector<double> >& cons_vars);
}

// velocity and pressure
vector< vector<double> > PrimitiveVariables::get_primitive_variables( vector< vector<double> >& cons_vars ){
    vector< vector<double> > prim_vars;

    int size = cons_vars[0].size();

    vector<double> velocity(size);
    vector<double> pressure(size);

    for(int i=0 ; i<size ; i++){
        velocity[i] = cons_vars[1][i] / cons_vars[0][i];
        pressure[i] = (CTS::GAMMA-1)*(cons_vars[2][i] - velocity[i]*cons_vars[1][i]*0.5);
    }

    prim_vars.push_back(velocity);
    prim_vars.push_back(pressure);

    return prim_vars;
}
