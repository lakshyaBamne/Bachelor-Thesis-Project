/*
    * @author Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 1-Dimensional Central Upwind Scheme
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Calculate CU Numerical Flux for a given set of Conserved variable vectors
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<cmath>

#include "Constants.h"
#include "Utility.h"
#include "ExtendCells.h"

using namespace std;

namespace CTS = Constants;
namespace UTL = Utility;
namespace EXC = ExtendCells;

namespace CUNumericalFlux{
    // main calculator
    vector<vector<double>> get_cu_flux(vector<vector<double>>& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time);

    // utility or intermediate functions
    vector<vector<double>> get_plr(vector<double>& vars, pair<string,string>& initial_conditions);
    vector<double> get_slx(vector<double>& vars, pair<string,string>& initial_conditions);

    // extra functions
    vector<vector<double>> get_plr(vector<double>& vars, pair<string,string>& initial_conditions, string var);
    vector<double> get_slx(vector<double>& vars, pair<string,string>& initial_conditions, string var);

    // functions to calculate the PLR for primitive variables used while calculation of the Local Speeds of Propagation
    vector<vector<double>> primitive_plr(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr);
    vector<vector<double>> primitive_plr(vector<vector<double>>& u_plr ,vector<vector<double>>& rho_plr, vector<vector<double>>& E_plr);

    // function to get the CU Local Speeds of Propagation given the Piecewise Linear Reconstructions
    vector<vector<double>> get_cu_lsp(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);

    // function to get the flux vectors given the Conserved variables
    vector<vector<double>> get_flux(int num, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);

    double new_dt(vector<vector<double>>& lsp, double dx, double t,  pair<double,double>& time);

    // Functions used to calculate the Anti-diffusion term
    vector<vector<double>> CalculateStar(vector<vector<double>>& lsp, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);
    vector<vector<double>> CalculateADT(vector<vector<double>>& Star, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);

    // Function to get the CU Numerical Flux vector
    vector<vector<double>> get_cu_flux(vector<vector<double>>& lsp, vector<vector<double>>& ADTerm, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3);
}

// Function to get the CU Numerical flux given any set of Conserved variables
/**
*   @arg cons_vars : 2D vector containing the Conserved variables (rho, m, E)
*   @arg initial_conditions : Riemann Problem , Boundary Conditions
*   @arg dt : time step for next iteration
*   @arg dx : volume(here length) of one cell in the computational domain
*   @arg t : current time (updated after every iteration using SSPRK)
**/
vector<vector<double>> CUNumericalFlux::get_cu_flux( vector<vector<double>>& cons_vars, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time){
    // STEP-1 Get the Piecewise Linear Reconstruction for the conserved variables
    // vector<vector<double>> rho_plr = get_plr(cons_vars[0], initial_conditions);
    // vector<vector<double>> m_plr = get_plr(cons_vars[1], initial_conditions);    
    // vector<vector<double>> E_plr = get_plr(cons_vars[2], initial_conditions);

    // new code to accomodate more boundary conditions
    vector< vector<double> > rho_plr = get_plr( cons_vars[0], initial_conditions, "density" );
    vector< vector<double> > m_plr = get_plr( cons_vars[1], initial_conditions, "momentum" );
    vector< vector<double> > E_plr = get_plr( cons_vars[2], initial_conditions, "energy" );

    // Now we need to calculate the One-Sided Local speeds of propagation (for CU schemes)
    vector<vector<double>> lsp = get_cu_lsp(rho_plr, m_plr, E_plr);
    
    dt = new_dt(lsp, dx, t, time); // updating time step for the next iteration

    // Now we need to get the Anti-Diffusion term 
    // -> for that first we need the term U* from equation 2.4 of the Research paper
    // -> to calculate U* we need the Flux vectors which can be calculated from the Conserved variables
    vector<vector<double>> f1 = get_flux(1, rho_plr, m_plr, E_plr);
    vector<vector<double>> f2 = get_flux(2, rho_plr, m_plr, E_plr);
    vector<vector<double>> f3 = get_flux(3, rho_plr, m_plr, E_plr);

    // Anti-diffusion vector (used only when required)
    vector<vector<double>> star = CalculateStar(lsp, f1, f2, f3, rho_plr, m_plr, E_plr);
    
    vector<vector<double>> ADTerm = CalculateADT(star, rho_plr, m_plr, E_plr);

    // Now we can calulate the CU Numerical Flux vector using which the next iteration is obtained
    vector<vector<double>> cu_flux = get_cu_flux(lsp, ADTerm, rho_plr, m_plr, E_plr, f1, f2, f3);

    return cu_flux;
}

// Function to calcualte dt for the next iteration using CFL conditions
double CUNumericalFlux::new_dt(vector<vector<double>>& lsp, double dx, double t,  pair<double,double>& time){
    double amax=0.0;

    int size = lsp[0].size();

    for(int i=0 ; i<size ; i++){
        amax = max( amax , max( lsp[0][i], -1*lsp[1][i] ) );
    }

    double newdt = (CTS::CFL * dx)/amax;

    if( t + newdt > time.second ){
        newdt = time.second - t;
    }

    return newdt;
}

// Function to calculate the Piecewise Linear Reconstruction of Conserved variables
vector<vector<double>> CUNumericalFlux::get_plr(vector<double>& vars, pair<string,string>& initial_conditions){
    // Calculate the slopes using the min mod limiter
    vector<double> slx = get_slx(vars, initial_conditions);
    
    vector<vector<double>> plr;

    vector<double> plr_East(vars.size()-2);
    vector<double> plr_West(vars.size()-2);

    for(int i=1 ; i<vars.size()-1 ; i++){
        plr_East[i-1] = vars[i] + slx[i]/2;
        plr_West[i-1] = vars[i] - slx[i]/2;
    }

    // extend cells
    EXC::extend_cells_ind(initial_conditions.second, plr_East);
    EXC::extend_cells_ind(initial_conditions.second, plr_West);

    plr.push_back(plr_East);
    plr.push_back(plr_West);

    return plr;
}

vector<double> CUNumericalFlux::get_slx(vector<double>& vars, pair<string,string>& initial_conditions){
    vector<double> slx(vars.size()-2);

    for(int i=1 ; i<vars.size()-1 ; i++){
        double val1 = CTS::THETA * ( vars[i+1]-vars[i] );
        double val2 = ( vars[i+1]-vars[i-1] )/2;
        double val3 = CTS::THETA * ( vars[i]-vars[i-1] );

        slx[i-1] = UTL::minmod( val1, val2, val3 );
    }

    // extend cell
    EXC::extend_cells_ind(initial_conditions.second , slx);

    return slx;
}

// Function to find the Piecewise Linear Reconstruction for velocity
vector<vector<double>> CUNumericalFlux::primitive_plr(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr){
    vector<vector<double>> u_plr;
    
    int size = rho_plr[0].size();

    vector<double> u_E(size);
    vector<double> u_W(size);

    for(int i=0 ; i<size ; i++){
        u_E[i] = m_plr[0][i] / rho_plr[0][i];
        u_W[i] = m_plr[1][i] / rho_plr[1][i];
    }

    // Ghost values are not required for reconstruction of Primitive variables

    u_plr.push_back(u_E);
    u_plr.push_back(u_W);

    return u_plr;
}

// Function to find the Piecewise Linear Reconstruction for Pressure
vector<vector<double>> CUNumericalFlux::primitive_plr(vector<vector<double>>& u_plr, vector<vector<double>>& rho_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> p_plr;

    int size = u_plr[0].size();

    vector<double> p_E(size);
    vector<double> p_W(size);

    for(int i=0 ; i<size ; i++){
        p_E[i] = (CTS::GAMMA-1)*( E_plr[0][i] - ( (rho_plr[0][i]*u_plr[0][i]*u_plr[0][i] )/2 ) );
        p_W[i] = (CTS::GAMMA-1)*( E_plr[1][i] - ( (rho_plr[1][i]*u_plr[1][i]*u_plr[1][i] )/2 ) );
    }

    // Ghost values are not required for reconstruction of Primitive variables

    p_plr.push_back(p_E);
    p_plr.push_back(p_W);

    return p_plr;
}

// CU Local speeds of propagation a+ and a-
//! a+ and a- are calculated for all the cells other than the right ghost value cell
vector<vector<double>> CUNumericalFlux::get_cu_lsp(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> lsp;

    // To apply the formula we need u,p,rho which need to be calculated using vector cv
    vector<vector<double>> u_plr = primitive_plr(rho_plr, m_plr);
    vector<vector<double>> p_plr = primitive_plr(u_plr, rho_plr, E_plr);

    int size = rho_plr[0].size()-1;

    vector<double> ap(size); // a+
    vector<double> am(size); // a-

    for(int i=0 ; i<size ; i++){
        double speed_sound_E =  sqrtf( (CTS::GAMMA*p_plr[0][i])/rho_plr[0][i] );
        double speed_sound_W = sqrtf( (CTS::GAMMA*p_plr[1][i+1])/rho_plr[1][i+1] );

        // cout << "speed of sound (E) : " << speed_sound_E << endl;
        // cout << "speed of sound (W) : " << speed_sound_W << endl;

        //! check if this has a problem
        // a+
        ap[i] = max( 0.0 , max( (u_plr[0][i] + speed_sound_E) , (u_plr[1][i+1] + speed_sound_W) ) );

        // a-
        am[i] = min( 0.0 , min( (u_plr[0][i] - speed_sound_E) , (u_plr[1][i+1] - speed_sound_W) ) );
    }

    // output
    lsp.push_back(ap);
    lsp.push_back(am);

    return lsp;
}

// Get the flux vectors given the Conserved variables
vector<vector<double>> CUNumericalFlux::get_flux(int num, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> flux;
    
    vector<vector<double>> u_plr = primitive_plr(rho_plr, m_plr);
    vector<vector<double>> p_plr = primitive_plr(u_plr, rho_plr, E_plr);

    int size = u_plr[0].size();

    vector<double> fE(size);    
    vector<double> fW(size);    

    switch (num){
        //! f1 = rho*u
        case 1:
            for(int i=0 ; i<size ; i++){
                fE[i] = m_plr[0][i];
                fW[i] = m_plr[1][i];
            }

            break;
        //! f2 = m*u + p
        case 2:
            for(int i=0 ; i<size ; i++){
                fE[i] = m_plr[0][i]*u_plr[0][i] + p_plr[0][i];
                fW[i] = m_plr[1][i]*u_plr[1][i] + p_plr[1][i];
            }

            break;
        //! f3 = u*(E + p)
        case 3:
            for(int i=0 ; i<size ; i++){
                fE[i] = u_plr[0][i]*( E_plr[0][i] + p_plr[0][i] );
                fW[i] = u_plr[1][i]*( E_plr[1][i] + p_plr[1][i] );
            }

            break;
    }

    flux.push_back(fE);
    flux.push_back(fW);

    return flux;
}

// Function to Calculate the U* values used for the Anti-diffusion term
vector<vector<double>> CUNumericalFlux::CalculateStar(vector<vector<double>>& lsp, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr ){
    vector<vector<double>> Star;

    int size = lsp[0].size();

    vector<double> rho_star(size);
    vector<double> m_star(size);
    vector<double> E_star(size);

    for(int i=0 ; i<size ; i++){
        rho_star[i] = ( lsp[0][i]*rho_plr[1][i+1] - lsp[1][i]*rho_plr[0][i] - (f1[1][i+1]-f1[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
        m_star[i] = ( lsp[0][i]*m_plr[1][i+1] - lsp[1][i]*m_plr[0][i] - (f2[1][i+1]-f2[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
        E_star[i] = ( lsp[0][i]*E_plr[1][i+1] - lsp[1][i]*E_plr[0][i] - (f3[1][i+1]-f3[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
    }

    Star.push_back(rho_star);
    Star.push_back(m_star);
    Star.push_back(E_star);

    return Star;
}

// Function to calculate the Anti-diffusion term
vector<vector<double>> CUNumericalFlux::CalculateADT(vector<vector<double>>& Star, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    int size = Star[0].size();

    vector<vector<double>> ADTerm;

    vector<double> d_rho(size);
    vector<double> d_m(size);
    vector<double> d_E(size);

    for(int i=0 ; i<size ; i++){
        d_rho[i] = UTL::minmod( rho_plr[1][i+1]-Star[0][i] , Star[0][i] - rho_plr[0][i] );
        d_m[i] = UTL::minmod( m_plr[1][i+1]-Star[1][i] , Star[1][i] - m_plr[0][i] );
        d_E[i] = UTL::minmod( E_plr[1][i+1]-Star[2][i] , Star[2][i] - E_plr[0][i] );
    }

    ADTerm.push_back(d_rho);
    ADTerm.push_back(d_m);
    ADTerm.push_back(d_E);

    return ADTerm;
}

vector<vector<double>> CUNumericalFlux::get_cu_flux(vector<vector<double>>& lsp, vector<vector<double>>& ADTerm, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3){
    // value of the flux depends on the difference between the Local Speeds of Propagation
    vector<vector<double>> cu_flux;

    int size = lsp[0].size();

    vector<double> cu_f1(size);
    vector<double> cu_f2(size);
    vector<double> cu_f3(size);

    for(int i=0 ; i<size ; i++){
        if( lsp[0][i]-lsp[1][i] > CTS::EPSILON ){
            // a+ - a- appears in the denominator so it cant be too small

            cu_f1[i] = ( lsp[0][i]*f1[0][i] - lsp[1][i]*f1[1][i+1] + lsp[0][i]*lsp[1][i]*( rho_plr[1][i+1] - rho_plr[0][i] - ADTerm[0][i] ) ) / ( lsp[0][i]-lsp[1][i] );
            cu_f2[i] = ( lsp[0][i]*f2[0][i] - lsp[1][i]*f2[1][i+1] + lsp[0][i]*lsp[1][i]*( m_plr[1][i+1] - m_plr[0][i] - ADTerm[1][i] ) ) / ( lsp[0][i]-lsp[1][i] );
            cu_f3[i] = ( lsp[0][i]*f3[0][i] - lsp[1][i]*f3[1][i+1] + lsp[0][i]*lsp[1][i]*( E_plr[1][i+1] - E_plr[0][i] - ADTerm[2][i] ) ) / ( lsp[0][i]-lsp[1][i] );
        }
        else{
            cu_f1[i] = 0.5 * ( f1[0][i] + f1[1][i] );
            cu_f2[i] = 0.5 * ( f2[0][i] + f2[1][i] );
            cu_f3[i] = 0.5 * ( f3[0][i] + f3[1][i] );
        }
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);

    return cu_flux;
}




/*
    Extra functions
*/
vector<vector<double>> CUNumericalFlux::get_plr(vector<double>& vars, pair<string,string>& initial_conditions, string var){
    // Calculate the slopes using the min mod limiter
    vector<double> slx = get_slx(vars, initial_conditions, var);
    
    vector<vector<double>> plr;

    vector<double> plr_East(vars.size()-2);
    vector<double> plr_West(vars.size()-2);

    for(int i=1 ; i<vars.size()-1 ; i++){
        plr_East[i-1] = vars[i] + slx[i]/2;
        plr_West[i-1] = vars[i] - slx[i]/2;
    }

    // extend cells
    EXC::extend_cells_ind(initial_conditions.second, plr_East, var);
    EXC::extend_cells_ind(initial_conditions.second, plr_West, var);

    plr.push_back(plr_East);
    plr.push_back(plr_West);

    return plr;
}

vector<double> CUNumericalFlux::get_slx(vector<double>& vars, pair<string,string>& initial_conditions, string var){
    vector<double> slx(vars.size()-2);

    for(int i=1 ; i<vars.size()-1 ; i++){
        double val1 = CTS::THETA * ( vars[i+1]-vars[i] );
        double val2 = ( vars[i+1]-vars[i-1] )/2;
        double val3 = CTS::THETA * ( vars[i]-vars[i-1] );

        slx[i-1] = UTL::minmod( val1, val2, val3 );
    }

    // extend cell
    EXC::extend_cells_ind(initial_conditions.second , slx, var);

    return slx;
}


