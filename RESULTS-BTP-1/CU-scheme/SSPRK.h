/*
    * @author Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 1-Dimensional Central Upwind Scheme
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Namespace to calculate the Next iteration values for the Conserved vectors
    -> given the current conserved vectors and the CU Numerical Flux so calculated
    -> 3-stage SSPRK scheme is used to do this calculation
*/
#pragma once

#include<iostream>
#include<vector>
#include<utility>

#include "ExtendCells.h"
#include "CUNumericalFlux.h"

using namespace std;

namespace CUF = CUNumericalFlux;
namespace EXC = ExtendCells;

namespace SSPRK{
    vector<vector<double>> get_next_cons_vars(vector<vector<double>>& cons_vars, vector<vector<double>>& cu_flux, double LMBDA, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time);
}

vector<vector<double>> SSPRK::get_next_cons_vars(vector<vector<double>>& cons_vars, vector<vector<double>>& cu_flux, double LAMBDA, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time){
    // we have to implement the 3-stage SSP-RK scheme to find new values of Conserved variables
    // given the old values and CU Numerical Flux calculated using the algorithm provided in the paper

    const double original_dt = dt; // we dont want to change the value of dt inside this function

    int size = cons_vars[0].size()-2;

    // 1) Calculate U1
    vector<vector<double>> U1;

    vector<double> U11(size);
    vector<double> U12(size);
    vector<double> U13(size);

    for(int i=1 ; i<=size ; i++){
        U11[i-1] = cons_vars[0][i] - LAMBDA * (cu_flux[0][i] - cu_flux[0][i-1]);
        U12[i-1] = cons_vars[1][i] - LAMBDA * (cu_flux[1][i] - cu_flux[1][i-1]);
        U13[i-1] = cons_vars[2][i] - LAMBDA * (cu_flux[2][i] - cu_flux[2][i-1]);
    }

    U1.push_back(U11);
    U1.push_back(U12);
    U1.push_back(U13);

    EXC::extend_cells(initial_conditions.second, U1);

    // calculating CU Numerical Flux for U1
    // ! we don't want to calculate a new time step here so we have to keep a copy of the old one
    vector<vector<double>> cu_flux_u1 = CUF::get_cu_flux(U1, initial_conditions, dt, dx, t, time);
    dt = original_dt;

    // 2) Calculate U2
    vector<vector<double>> U2;

    vector<double> U21(size);
    vector<double> U22(size);
    vector<double> U23(size);

    for(int i=1 ; i<=size ; i++){
        U21[i-1] = ( 3*cons_vars[0][i] + U1[0][i] - LAMBDA*(cu_flux_u1[0][i] - cu_flux_u1[0][i-1]) ) / 4;
        U22[i-1] = ( 3*cons_vars[1][i] + U1[1][i] - LAMBDA*(cu_flux_u1[1][i] - cu_flux_u1[1][i-1]) ) / 4;
        U23[i-1] = ( 3*cons_vars[2][i] + U1[2][i] - LAMBDA*(cu_flux_u1[2][i] - cu_flux_u1[2][i-1]) ) / 4;
    }

    U2.push_back(U21);
    U2.push_back(U22);
    U2.push_back(U23);

    EXC::extend_cells(initial_conditions.second, U2);

    // calculating CU Numerical Flux for U2
    // ! we don't want to calculate a new time step here so we have to keep a copy of the old one
    vector<vector<double>> cu_flux_u2 = CUF::get_cu_flux(U2, initial_conditions, dt, dx, t, time);
    dt = original_dt;

    // 3) Calculate U(n+1)
    vector<vector<double>> Unew;

    vector<double> Un1(size);
    vector<double> Un2(size);
    vector<double> Un3(size);

    for(int i=1 ; i<=size ; i++){
        Un1[i-1] = ( cons_vars[0][i] + 2*U2[0][i] - 2*LAMBDA*( cu_flux_u2[0][i] - cu_flux_u2[0][i-1] ) ) / 3;
        Un2[i-1] = ( cons_vars[1][i] + 2*U2[1][i] - 2*LAMBDA*( cu_flux_u2[1][i] - cu_flux_u2[1][i-1] ) ) / 3;
        Un3[i-1] = ( cons_vars[2][i] + 2*U2[2][i] - 2*LAMBDA*( cu_flux_u2[2][i] - cu_flux_u2[2][i-1] ) ) / 3;
    }

    Unew.push_back(Un1);
    Unew.push_back(Un2);
    Unew.push_back(Un3);

    EXC::extend_cells(initial_conditions.second, Unew);

    return Unew;

}


