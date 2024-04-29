/*
    Code for 1st Order Local Lax Friedrichs Scheme
*/

#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<fstream>
#include<sstream>
#include<math.h>
#include<cmath>
#include<chrono>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define pdd pair<double,double>
#define pss pair<string,string>

namespace CTS{
    //! Global constants
    double GAMMA = 1.4;
    double THETA = 1.3;
    double EPS = 1.0E-12;

    //! CFL number
    double CFL = 0.9; // 1st order schemes

    //! domain endpoints
    pdd domainX = {0, 1};
    // pdd domainX = {-5, 5}; // LAX

    //! time endpoints
    // pdd time = {0, 2}; // MCW
    // pdd time = {0, 0.038}; // BLW
    // pdd time = {0, 1.3}; // LAX
    pdd time = {0, 0.2}; // TORO-1
    // pdd time = {0, 0.012}; // TORO-3
    // pdd time = {0, 0.035}; // TORO-4

    //! problem 
    // string PROBLEM = "MCW";
    // string PROBLEM = "BLW";
    // string PROBLEM = "LAX";
    string PROBLEM = "TORO-1";
    // string PROBLEM = "TORO-2";
    // string PROBLEM = "TORO-3";

    string BC = "FREE";
    // string BC = "REFLECTIVE";

    //! number of finite volume cells
    // int N = 100;
    // int N = 1000;
    // int N = 2000;
    // int N = 3000;
    // int N = 4000;
    int N = 5000;
    // int N = 10000;
    // int N = 20000;

    int Nref = 4000;

    // length of an individual finite volume cell
    double dx = (domainX.second-domainX.first)/N;
}

#define vd vector<double>
#define vvd vector< vector<double> >
#define pdd pair<double,double>
#define pss pair<string,string>

// Function definitions
vd make_grid(pdd domainX, int N);
vvd initialize_conserved_variables(string PROBLEM, string BC, vd& grid, vd& rho, vd& u, vd& p, vd& a);
void extend_cells(string BC, vvd& U, vd& rho, vd& u, vd& p, vd& a);
void write_grid(vd& grid, string file);
void write_density(vd& density, string file);

void update_conserved_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a, double dx, double& dt);
vvd get_flux_vectors(vvd& U, vd& rho, vd& u, vd& p, vd& a);
vvd get_lf_flux(vvd& U, vvd& F, vd& u, vd& a, double dt);
void run_lf_scheme(string MODE);
double minmod(double a, double b, double c);
void update_primitive_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a);

int main(){
    run_lf_scheme("NORMAL");
    // run_lf_scheme("REF");

    return 0;
}

//! Function bodies

// Function to initialize a computational grid
vd make_grid(pdd domainX, int N){
    double dx = (domainX.second-domainX.first)/N;

    vd grid(N+2);

    for(int i=1 ; i<=N ; i++){
        grid[i] = domainX.first + (i - 0.5)*dx;
    }

    return grid;
}

// Function to initialize the conserved variables
vvd initialize_conserved_variables(string PROBLEM, string BC, vd& grid, vd& rho, vd& u, vd& p, vd& a){
    
    int N = rho.size()-2;

    if(PROBLEM == "MCW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.3){
                rho[i] = 1.4;
                u[i] = 0.1;
                p[i] = 1;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0.1;
                p[i] = 1;
            }
        }
    }
    else if(PROBLEM == "BLW"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.1){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 1000;
            }
            else if(grid[i]>=0.1 && grid[i]<=0.9){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 0.01;
            }
            else{
                rho[i] = 1;
                u[i] = 0;
                p[i] = 100;
            }
        }
    }
    else if(PROBLEM == "LAX"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0){
                rho[i] = 0.445;
                u[i] = 0.698;
                p[i] = 3.528;
            }
            else{
                rho[i] = 0.500;
                u[i] = 0;
                p[i] = 0.571;
            }
        }
    }
    else if(PROBLEM == "TORO-1"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.3){
                rho[i] = 1.0;
                u[i] = 0.75;
                p[i] = 1.0;
            }
            else{
                rho[i] = 0.125;
                u[i] = 0;
                p[i] = 0.1;
            }
        }
    }
    else if(PROBLEM == "TORO-2"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.5){
                rho[i] = 1.0;
                u[i] = 0;
                p[i] = 1000;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0;
                p[i] = 0.01;
            }
        }
    }
    else if(PROBLEM == "TORO-3"){
        for(int i=1 ; i<=N ; i++){
            if(grid[i] < 0.4){
                rho[i] = 5.99924;
                u[i] = 19.5975;
                p[i] = 460.894;
            }
            else{
                rho[i] = 5.99242;
                u[i] = -6.19633;
                p[i] = 46.0950;
            }
        }
    }
    else{
        cout << "---LOG--- Please enter correct problem---" << endl;
    }

    // now initialize the conserved variables
    vvd U(3, vd(N+2, 0));

    for(int i=1 ; i<=N ; i++){
        U[0][i] = rho[i];
        U[1][i] = rho[i]*u[i];
        U[2][i] = (p[i]/(CTS::GAMMA-1)) + (0.5*rho[i]*pow(u[i],2));

        a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);
    }

    return U;
}

// Function to extend cells for the conserved variables
void extend_cells(string BC, vvd& U, vd& rho, vd& u, vd& p, vd& a){
    int N = rho.size()-2;

    if(BC == "FREE"){
        // Density
        rho[0] = rho[1];
        rho[N+1] = rho[N];

        // Velocity
        u[0] = u[1];
        u[N+1] = u[N];

        // Presure 
        p[0] = p[1];
        p[N+1] = p[N];

        // Acoustic speed
        a[0] = a[1];
        a[N+1] = a[N];

        // Conserved variables
        for(int i=0 ; i<3 ; i++){
            U[i][0] = U[i][1];
            U[i][N+1] = U[i][N];
        }

    }
    else if(BC == "REFLECTIVE"){
        // Density
        rho[0] = rho[1];
        rho[N+1] = rho[N];

        // Velocity
        u[0] = -u[1];
        u[N+1] = -u[N];

        // Pressure
        p[0] = p[1];
        p[N+1] = p[N];

        // Acoustic speed
        a[0] = a[1];
        a[N+1] = a[N];

        // Conserved variable-1
        U[0][0] = U[0][1];
        U[0][N+1] = U[0][N];

        // Conserved variable-2
        U[1][0] = -U[1][1];
        U[1][N+1] = -U[1][N];

        // Conserved variable-3
        U[2][0] = U[2][1];
        U[2][N+1] = U[2][N];

    }
    else{
        cout << "---LOG--- Please enter correct boundary conditions---" << endl;
    }
}

// Function to write the computational grid to external file
void write_grid(vd& grid, string file){
    stringstream ss;
    for(int i=1 ; i<grid.size()-1 ; i++){
        ss << grid[i] << " ";
    }
    
    ofstream fout;
    fout.open(file, ios::app);

    fout << ss.str() << endl;

    fout.close();
}

// Function to write the density to an external file
void write_density(vd& density, string file){
    stringstream ss;
    for(int i=1 ; i<density.size()-1 ; i++){
        ss << density[i] << " ";
    }

    ofstream fout;
    fout.open(file, ios::app);

    fout << ss.str() << endl;

    fout.close();
}

// Function to update the conserved variables using the LF Scheme
void update_conserved_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a, double dx, double& dt){
    int N = U[0].size()-2;

    // first we need to calculate the time step dt
    double amax;

    for(int i=1 ; i<=N ; i++){
        amax = max(
            amax,
            abs(u[i])+a[i]
        );
    }

    //! update time step
    dt = CTS::CFL*dx/amax;

    // first we need to extend cells for the variables to calculate intercell fluxes
    extend_cells(CTS::BC, U, rho, u, p, a);    

    //! Calculate the Fluxes at the finite volume grid points
    vvd F = get_flux_vectors(U, rho, u, p, a);

    //! Calculate the Lax-Friedrichs Fluxes
    vvd LF = get_lf_flux(U, F, u, a, dt);

    //! Update the variables using the calculated fluxes using Euler Forward Differences
    for(int i=1 ; i<=CTS::N ; i++){
        U[0][i] = U[0][i] - (dt/CTS::dx)*(LF[0][i]-LF[0][i-1]);
        U[1][i] = U[1][i] - (dt/CTS::dx)*(LF[1][i]-LF[1][i-1]);
        U[2][i] = U[2][i] - (dt/CTS::dx)*(LF[2][i]-LF[2][i-1]);
    }

    // now we can update the other variables 
    for(int i=1 ; i<=CTS::N ; i++){
        rho[i] = U[0][i];
        u[i] = U[1][i]/U[0][i];
        p[i] = (CTS::GAMMA-1)*(U[2][i] - rho[i]*pow(u[i],2));
        a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);
    }

}

// Function to calculate the Flux vectors at the grid points
vvd get_flux_vectors(vvd& U, vd& rho, vd& u, vd& p, vd& a){
    vvd F(3, vd(CTS::N+2));

    for(int i=0 ; i<=CTS::N+1 ; i++){
        F[0][i] = U[1][i];
        F[1][i] = U[0][i]*pow(u[i],2) + p[i];
        F[2][i] = u[i]*(U[2][i] + p[i]);
    }

    return F;
}

// Function to calculate the LF Fluxes
vvd get_lf_flux(vvd& U, vvd& F, vd& u, vd& a, double dt){
    vvd LF(3, vd(CTS::N+1));

    for(int i=0 ; i<=CTS::N ; i++){
        // for rusanov's scheme we need to find the intermediate parameter
        double alpha = max(
            abs(u[i])+a[i],
            abs(u[i+1])+a[i+1]
        );

        LF[0][i] = 0.5*(F[0][i]+F[0][i+1]) - 0.5*alpha*(U[0][i+1]-U[0][i]);
        LF[1][i] = 0.5*(F[1][i]+F[1][i+1]) - 0.5*alpha*(U[1][i+1]-U[1][i]);
        LF[2][i] = 0.5*(F[2][i]+F[2][i+1]) - 0.5*alpha*(U[2][i+1]-U[2][i]);
    }

    return LF;
}

void run_lf_scheme(string MODE){
    if(MODE == "NORMAL"){
        //! initialize computational grid
        vd grid = make_grid(CTS::domainX, CTS::N);

        //! output the computational grid to external file
        write_grid(grid, MODE+"grid.txt");

        //! initialize the conserved variables vector based on the given problem
        vd rho(CTS::N+2), u(CTS::N+2), p(CTS::N+2), a(CTS::N+2);
        vvd U = initialize_conserved_variables(CTS::PROBLEM, CTS::BC, grid, rho, u, p, a);

        double dt;
        double t=CTS::time.first;

        auto start = chrono::high_resolution_clock::now();

        write_density(U[0], MODE+"density.txt");

        cout << "Serial scheme running" << endl;

        while(t < CTS::time.second){
            // log message
            // cout << "t=" << t << "|dt=" << dt << endl;

            // write to output file

            update_conserved_variables(U, rho, u, p, a, CTS::dx, dt);
            t += dt;
        }

        write_density(U[0], MODE+"density.txt");

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = end - start;

        cout << "Scheme run time (serial): " << duration.count() << endl;
    }
    else if(MODE == "REF"){
        //! initialize computational grid
        vd grid = make_grid(CTS::domainX, CTS::Nref);

        //! output the computational grid to external file
        write_grid(grid, MODE+"grid.txt");

        //! initialize the conserved variables vector based on the given problem
        vd rho(CTS::N+2), u(CTS::N+2), p(CTS::N+2), a(CTS::N+2);
        vvd U = initialize_conserved_variables(CTS::PROBLEM, CTS::BC, grid, rho, u, p, a);

        double dt;
        double t=CTS::time.first;

        while(t < CTS::time.second){
            // log message
            cout << "t=" << t << "|dt=" << dt << endl;

            // write to output file
            write_density(U[0], MODE+"density.txt");

            update_conserved_variables(U, rho, u, p, a, CTS::dx, dt);
            t += dt;
        
        }
    }
    else{
        cout << "---ERROR--- Please select a correct run mode---" << endl;
    }
}

// Minmod limiter used to supress oscillations in slope reconstructions
double minmod(double a, double b, double c){
    if( min(a, min(b,c)) > 0 ){
        return min(a, min(b,c));
    }
    else if( max(a, max(b, c)) < 0 ){
        return max(a, max(b, c));
    }
    else{
        return 0;
    }
}

void update_primitive_variables(vvd& U, vd& rho, vd& u, vd& p, vd& a){

    for(int i=1 ; i<=CTS::N ; i++){
        rho[i] = U[0][i];
        u[i] = U[1][i]/U[0][i];
        p[i] = (CTS::GAMMA-1)*(U[2][i] - 0.5*U[0][i]*pow(u[i],2));
        a[i] = sqrt(CTS::GAMMA*p[i]/U[0][i]);
    }

}


