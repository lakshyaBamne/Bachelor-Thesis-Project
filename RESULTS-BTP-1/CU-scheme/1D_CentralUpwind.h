/*
    * @author Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 1-Dimensional Central Upwind Scheme
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/
#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<string>
#include<cmath>
#include<sstream>
#include<fstream>

#include "Constants.h"
#include "Utility.h"
#include "InitializeRiemannProblem.h"
#include "ExtendCells.h"
#include "CUNumericalFlux.h"
#include "PrimitiveVariables.h"
#include "SSPRK.h"
#include "OutputResult.h"

namespace IRP = InitRiemannProblem;
namespace UTL = Utility;
namespace EXC = ExtendCells;
namespace CTS = Constants;
namespace CUF = CUNumericalFlux;
namespace PRV = PrimitiveVariables;
namespace SRK = SSPRK;
namespace OPF = OutputResult;

#define ll long long int

using namespace std;

/*
    Functions to get input for starting the simulations
*/
class GetInput{

public: // variables
    pair<double,double> domain;
    pair<double,double> time;
    ll Nx;
    pair<string,string> initial_conditions;

    double dt;
    double t;
    double dx;

public: // methods
    GetInput(){
        // show to the user information about the lubrary and it's basic usage
        show_info();

        // take input for the starting state
        cout << "|---------------------------------------INPUT----------------------------------------|" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;

        cout << "Domain start : ";
        cin >> domain.first;

        cout << "Domain end : ";
        cin >> domain.second;

        cout << "Grid points in domain : ";
        cin >> Nx;

        cout << "Initial time : ";
        cin >> time.first;

        cout << "Final time : ";
        cin >> time.second;

        cout << "Enter the Riemann Problem : ";
        cin >> initial_conditions.first;

        cout << "Enter the Boundary Conditions : ";
        cin >> initial_conditions.second;

        cout << "|------------------------------------------------------------------------------------|" << endl;

        t = time.first;
        dx = ( domain.second - domain.first ) / Nx;

        stringstream ss1, ss2;
        ss1 << Nx;

        UTL::export_string(initial_conditions.first);
        UTL::export_string(ss1.str());

        ss2 << domain.first << " " << domain.second;
        UTL::export_string(ss2.str());
    }

    // function to show information about the problems solvable by the system
    void show_info(){
        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                       1-Dimensional Central Upwind Scheme                          |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                  Riemann problems solvable using this library                      |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;

        cout << "-> [MCW] Slowly moving isolated contact discontinuity" << endl;
        cout << "-> [SCW] Stationary Contact Wave and Travelling shock and rarefaction" << endl;
        cout << "-> [LAX] Lax Problem" << endl;
        cout << "-> [BLW] Blastwave Problem" << endl;
        cout << "-> [SDW] Shock-density wave interaction problem" << endl;
        cout << "-> [SEW] Shock-entropy wave interaction problem" << endl;

        cout << "+------------------------------------------------------------------------------------+" << endl;
        cout << "|                            Boundary Conditions supported                           |" << endl;
        cout << "+------------------------------------------------------------------------------------+" << endl;
        
        cout << "-> [FREE] Neumann Boundary Conditions" << endl;
        cout << "-> [REFLECTIVE] Reflective/Solidwall Boundary Conditions" << endl;
        cout << "-> [PERIODIC] Periodic Boundary Conditions" << endl;
        cout << "-> [RFREE] Reflect Free Boundary Conditions" << endl;

        cout << "+------------------------------------------------------------------------------------+" << endl;
    }

    /*
        Function to run complete iterations for the Central Upwind Numerical Scheme
        -> store data for all the iterations
    */
    void run_cu_scheme_complete(string result){

        vector<vector<double>> cons_vars = IRP::get_conserved_variables(initial_conditions, domain, Nx, result);

        string file_density = result + "/Density.txt";
        string file_momentum = result + "/Momentum.txt";
        string file_energy = result + "/Energy.txt";

        string file_velocity = result + "/Velocity.txt";
        string file_pressure = result + "/Pressure.txt";

        while( t < time.second ){
            // Log
            cout << "t = " << t << " | dt = " << dt << endl;

            //! Writing output to a file for plotting later
            OPF::write_vector(cons_vars[0], file_density);
            OPF::write_vector(cons_vars[1], file_momentum);
            OPF::write_vector(cons_vars[2], file_energy);

            // convert the conserved variables to primitive variables for plotting
            vector< vector<double> > prim_vars = PRV::get_primitive_variables(cons_vars);

            OPF::write_vector(prim_vars[0], file_velocity);
            OPF::write_vector(prim_vars[1], file_pressure);

            vector<vector<double>> cu_flux = CUF::get_cu_flux(cons_vars, initial_conditions, dt, dx, t, time);

            double LAMBDA = dt / dx;
            
            t = t+dt;

            vector<vector<double>> cons_vars_next = SRK::get_next_cons_vars(cons_vars, cu_flux, LAMBDA, initial_conditions, dt, dx, t, time);    
        
            // copy the new values in the old vector to be used in the next iteration
            for(int i=0 ; i<cons_vars[0].size() ; i++){
                cons_vars[0][i] = cons_vars_next[0][i];
                cons_vars[1][i] = cons_vars_next[1][i];
                cons_vars[2][i] = cons_vars_next[2][i];
            }
        }

    }

    void run_cu_scheme_partial(string result){
        vector<vector<double>> cons_vars = IRP::get_conserved_variables(initial_conditions, domain, Nx, result);
        
        string file_density = result + "/Density.txt";
        string file_momentum = result + "/Momentum.txt";
        string file_energy = result + "/Energy.txt";

        string file_velocity = result + "/Velocity.txt";
        string file_pressure = result + "/Pressure.txt";

        //! Writing output to a file for plotting later
        OPF::write_vector(cons_vars[0], file_density);
        OPF::write_vector(cons_vars[1], file_momentum);
        OPF::write_vector(cons_vars[2], file_energy);

        // convert the conserved variables to primitive variables for plotting
        vector< vector<double> > prim_vars = PRV::get_primitive_variables(cons_vars);

        OPF::write_vector(prim_vars[0], file_velocity);
        OPF::write_vector(prim_vars[1], file_pressure);

        while( t < time.second ){
            // Log
            cout << "t = " << t << " | dt = " << dt << endl;

            vector<vector<double>> cu_flux = CUF::get_cu_flux(cons_vars, initial_conditions, dt, dx, t, time);

            double LAMBDA = dt / dx;
            
            t = t+dt;

            vector<vector<double>> cons_vars_next = SRK::get_next_cons_vars(cons_vars, cu_flux, LAMBDA, initial_conditions, dt, dx, t, time);    
        
            // copy the new values in the old vector to be used in the next iteration
            for(int i=0 ; i<cons_vars[0].size() ; i++){
                cons_vars[0][i] = cons_vars_next[0][i];
                cons_vars[1][i] = cons_vars_next[1][i];
                cons_vars[2][i] = cons_vars_next[2][i];
            }
        }

        //! Writing output to a file for plotting later
        OPF::write_vector(cons_vars[0], file_density);
        OPF::write_vector(cons_vars[1], file_momentum);
        OPF::write_vector(cons_vars[2], file_energy);

        // convert the conserved variables to primitive variables for plotting
        prim_vars = PRV::get_primitive_variables(cons_vars);

        OPF::write_vector(prim_vars[0], file_velocity);
        OPF::write_vector(prim_vars[1], file_pressure);
    }


};

