#include<mpi.h>
#include<iostream>
#include<vector>
#include<string>
#include<utility>
#include<unordered_map>
#include<fstream>
#include<sstream>
#include<cmath>
#include<math.h>
#include<chrono>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>

#define MPI_DP MPI_DOUBLE_PRECISION
#define MPI_CW MPI_COMM_WORLD
#define MPI_SI MPI_STATUS_IGNORE

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
    // pdd time = {0, 0.012}; // TORO-2
    // pdd time = {0, 0.035}; // TORO-3

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

// function to write the output files to external text files
void write_grid(vd& grid, string file);
void write_density(vd& density, string file);

// parallel functions
void sync(int world_size, int rank);
void make_grid(int rank, vd& x, int single_size, int world_size);

// serial functions
double minmod(double a, double b, double c);
vvd initialize_conserved_variables(vd& grid, vd& rho, vd& u, vd& p, vd& a);
void extend_cells(vd& U1, vd& U2, vd& U3, vd& rho, vd& u, vd& p, vd& a);

int main(int argc, char* argv[]){
    // Initialize and get the command line arguments passed by the user
    MPI_Init(&argc, &argv);
    MPI_Status status;

    // get the rank of the current process that is running
    int rank, world_size;
    MPI_Comm_rank(MPI_CW, &rank);
    MPI_Comm_size(MPI_CW, &world_size);

    // variable declarations
    auto time_start = chrono::high_resolution_clock::now();
    auto time_end = chrono::high_resolution_clock::now();

    int sync_msg;
    int termination_msg = 4234532;
    int single_size = (int)(CTS::N/world_size);

    vd x(CTS::N), is;
    vd rho(CTS::N+2), u(CTS::N+2), p(CTS::N+2), a(CTS::N+2);
    vd U1(CTS::N+2), U2(CTS::N+2), U3(CTS::N+2);
    vd F1(CTS::N+2), F2(CTS::N+2), F3(CTS::N+2);

    vd lambda(CTS::N+1);
    vd LF_flux1(CTS::N+1), LF_flux2(CTS::N+1), LF_flux3(CTS::N+1);

    double amax=0, dt=0, t=CTS::time.first;

    // Initialize the computational grid (parallel)
    make_grid(rank, x, single_size, world_size);

    // sync processes
    sync(world_size, rank);

    if(rank == 0){
        cout << "---Parent running after initialization and sync---" << endl;
        write_grid(x, "parallel_grid.txt");

        // initialize conserved variables
        vvd U = initialize_conserved_variables(x, rho, u, p, a);
        for(int i=1 ; i<=CTS::N ; i++){
            U1[i] = U[0][i];
            U2[i] = U[1][i];
            U3[i] = U[2][i];
        }

        // extend the conserved variables
        extend_cells(U1, U2, U3, rho, u, p, a);

        for(int i=0 ; i<=CTS::N+1 ; i++){
            F1[i] = U2[i];
            F2[i] = (pow(U2[i],2)/U1[i]) + (CTS::GAMMA-1)*(U3[i] - 0.5*(pow(U2[i],2)/U1[i]));
            F3[i] = (U2[i]/U1[i])*( U3[i] + (CTS::GAMMA-1)*(U3[i] - 0.5*(pow(U2[i],2)/U1[i])) ); 
        }

        write_density(U1, "parallel_density.txt");
    }

    // sync processes
    sync(world_size, rank);

    auto start = chrono::high_resolution_clock::now();

    cout << "Starting Parallel Scheme" << endl;

    while(t < CTS::time.second){

        // Start the numerical scheme code 
        if(rank == 0){        
            // calculate the time step using the CFL conditions
            amax = 0;
            for(int i=1 ; i<=CTS::N ; i++){
                amax = max(
                    amax,
                    abs(u[i]) + a[i]
                );
            }

            dt = CTS::CFL*CTS::dx/amax;

            // cout << "t=" << t << " | dt=" << dt << endl;

            t += dt;
            
            if(t >= CTS::time.second){
                auto end = chrono::high_resolution_clock::now();
                chrono::duration<double> duration = end - start;

                cout << "Scheme run time (parallel): " << duration.count() << endl;

                write_density(U1, "parallel_density.txt");

                MPI_Abort(MPI_CW, 0);
            }

            // distribute equal work to all the nodes
            for(int i=1 ; i<world_size ; i++){
                MPI_Send(&U1[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
                MPI_Send(&U2[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
                MPI_Send(&U3[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);

                MPI_Send(&F1[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
                MPI_Send(&F2[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
                MPI_Send(&F3[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);

                MPI_Send(&u[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
                MPI_Send(&a[i*single_size], single_size+1, MPI_DP, i, i, MPI_CW);
            }

            // do the computation for the parent node
            for(int j=0 ; j<=single_size-1 ; j++){
                lambda[j] = max(
                    abs(u[j]) + a[j],
                    abs(u[j+1]) + a[j+1]
                );

                LF_flux1[j] = 0.5*(F1[j] + F1[j+1]) - 0.5*lambda[j]*(U1[j+1]-U1[j]);
                LF_flux2[j] = 0.5*(F2[j] + F2[j+1]) - 0.5*lambda[j]*(U2[j+1]-U2[j]);
                LF_flux3[j] = 0.5*(F3[j] + F3[j+1]) - 0.5*lambda[j]*(U3[j+1]-U3[j]);
            }

            // also do it for the last FVM cell
            lambda[CTS::N] = max(
                abs(u[CTS::N]) + a[CTS::N],
                abs(u[CTS::N+1]) + a[CTS::N+1]
            );

            LF_flux1[CTS::N] = 0.5*(F1[CTS::N] + F1[CTS::N+1]) - 0.5*lambda[CTS::N]*(U1[CTS::N+1]-U1[CTS::N]);
            LF_flux2[CTS::N] = 0.5*(F2[CTS::N] + F2[CTS::N+1]) - 0.5*lambda[CTS::N]*(U2[CTS::N+1]-U2[CTS::N]);
            LF_flux3[CTS::N] = 0.5*(F3[CTS::N] + F3[CTS::N+1]) - 0.5*lambda[CTS::N]*(U3[CTS::N+1]-U3[CTS::N]);

            // after the nodes do the computation gather data from them and update the results
            for(int i=1 ; i<world_size ; i++){
                MPI_Recv(&LF_flux1[i*single_size], single_size, MPI_DP, i, i, MPI_CW, MPI_SI);
                MPI_Recv(&LF_flux2[i*single_size], single_size, MPI_DP, i, i, MPI_CW, MPI_SI);
                MPI_Recv(&LF_flux3[i*single_size], single_size, MPI_DP, i, i, MPI_CW, MPI_SI);
            }

            // SERIAL CODE : To be parallelized

            // Now update the conserved variables
            for(int i=1 ; i<=CTS::N ; i++){
                U1[i] = U1[i] - (dt/CTS::dx)*(LF_flux1[i] - LF_flux1[i-1]);
                U2[i] = U2[i] - (dt/CTS::dx)*(LF_flux2[i] - LF_flux2[i-1]);
                U3[i] = U3[i] - (dt/CTS::dx)*(LF_flux3[i] - LF_flux3[i-1]);
            }

            extend_cells(U1, U2, U3, rho, u, p, a);

            // now we should also update the other variables for the next iteration
            // we need to update the flux variables for the parent process
            for(int i=0 ; i<=CTS::N+1 ; i++){
                F1[i] = U2[i];
                F2[i] = (pow(U2[i],2)/U1[i]) + (CTS::GAMMA-1)*(U3[i] - 0.5*(pow(U2[i],2)/U1[i]));
                F3[i] = (U2[i]/U1[i])*( U3[i] + (CTS::GAMMA-1)*(U3[i] - 0.5*(pow(U2[i],2)/U1[i])) );
            
                rho[i] = U1[i];
                u[i] = U2[i]/U1[i];
                p[i] = (CTS::GAMMA-1)*(U3[i] - 0.5*(pow(U2[i],2)/U1[i]));

                a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);    
            }


        }
        else{
            // receive data from the parent node
            MPI_Recv(&U1[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
            MPI_Recv(&U2[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
            MPI_Recv(&U3[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);

            MPI_Recv(&F1[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
            MPI_Recv(&F2[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
            MPI_Recv(&F3[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
        
            MPI_Recv(&u[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);
            MPI_Recv(&a[rank*single_size], single_size+1, MPI_DP, 0, rank, MPI_CW, MPI_SI);

            // make the required computations
            for(int j=rank*single_size ; j<=((rank+1)*single_size)-1 ; j++){
                lambda[j] = max(
                    abs(u[j]) + a[j],
                    abs(u[j+1]) + a[j+1]
                );

                LF_flux1[j] = 0.5*(F1[j] + F1[j+1]) - 0.5*lambda[j]*(U1[j+1]-U1[j]);
                LF_flux2[j] = 0.5*(F2[j] + F2[j+1]) - 0.5*lambda[j]*(U2[j+1]-U2[j]);
                LF_flux3[j] = 0.5*(F3[j] + F3[j+1]) - 0.5*lambda[j]*(U3[j+1]-U3[j]);
            }

            // return the results to the parent node (lambda and LF Fluxes)
            MPI_Send(&LF_flux1[rank*single_size], single_size, MPI_DP, 0, rank, MPI_CW);
            MPI_Send(&LF_flux2[rank*single_size], single_size, MPI_DP, 0, rank, MPI_CW);
            MPI_Send(&LF_flux3[rank*single_size], single_size, MPI_DP, 0, rank, MPI_CW);

        }
    }

    MPI_Finalize();
    return 0;
}

// function to sync all the processes at some point in the program
void sync(int world_size, int rank){
    int sync_msg;

    if(rank == 0){
        // send a message to all nodes then wait for a response from all nodes
        for(int i=1 ; i<world_size ; i++){
            MPI_Send(&sync_msg, 1, MPI_INT, i, 0, MPI_CW);
        }

        for(int i=1 ; i<world_size ; i++){
            MPI_Recv(&sync_msg, 1, MPI_INT, i, MPI_ANY_TAG, MPI_CW, MPI_SI);
        }
    }
    else{
        MPI_Recv(&sync_msg, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_CW, MPI_SI);
        MPI_Send(&sync_msg, 1, MPI_INT, 0, 0, MPI_CW);
    }
}

// Function to initialize the computational grid
void make_grid(int rank, vd& x, int single_size, int world_size){
    // Initialization for computational grid
    for(int j=rank*single_size ; j<=(((rank+1)*single_size) - 1) ; j++){
        x[j] = CTS::domainX.first + (j+0.5)*CTS::dx;
    }

    // now it sends the information to the parent
    if(rank == 0){
        for(int i=1 ; i<world_size ; i++){
            MPI_Recv(&x[i*single_size], single_size, MPI_DP, i, MPI_ANY_TAG, MPI_CW, MPI_SI);
        }
    }
    else{
        MPI_Send(&x[rank*single_size], single_size, MPI_DP, 0, 0, MPI_CW);
    }
}

void update_vars(){

}

// Function to write the computational grid to external file
void write_grid(vd& grid, string file){
    stringstream ss;
    for(int i=0 ; i<CTS::N ; i++){
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
    for(int i=1 ; i<=CTS::N ; i++){
        ss << density[i] << " ";
    }

    ofstream fout;
    fout.open(file, ios::app);

    fout << ss.str() << endl;

    fout.close();
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

vvd initialize_conserved_variables(vd& grid, vd& rho, vd& u, vd& p, vd& a){
    vvd U(3, vd(CTS::N+2));
    
    if(CTS::PROBLEM == "MCW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.3){
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
    else if(CTS::PROBLEM == "SCW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.8){
                rho[i] = 1;
                u[i] = -19.59745;
                p[i] = 1000;
            }
            else{
                rho[i] = 1;
                u[i] = -19.59745;
                p[i] = 0.01;
            }
        }
    }
    else if(CTS::PROBLEM == "BLW"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.1){
                rho[i] = 1;
                u[i] = 0;
                p[i] = 1000;
            }
            else if(grid[i-1]>=0.1 && grid[i-1]<=0.9){
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
    else if(CTS::PROBLEM == "LAX"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0){
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
    else if(CTS::PROBLEM == "TORO-1"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.3){
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
    else if(CTS::PROBLEM == "TORO-2"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.5){
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
    else if(CTS::PROBLEM == "TORO-3"){
        for(int i=1 ; i<=CTS::N ; i++){
            if(grid[i-1] < 0.4){
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
        cout << "---ERROR--- Please enter correct problem ---" << endl;
    }

    // initialize conserved variables
    for(int i=1 ; i<=CTS::N ; i++){
        U[0][i] = rho[i];
        U[1][i] = rho[i]*u[i];
        U[2][i] = ( p[i]/(CTS::GAMMA-1) ) + 0.5*rho[i]*pow(u[i],2);
        
        a[i] = sqrt(CTS::GAMMA*p[i]/rho[i]);
    }

    return U;
}

// Function to extend cells for the conserved variables
void extend_cells(vd& U1, vd& U2, vd& U3, vd& rho, vd& u, vd& p, vd& a){
    int N = CTS::N;

    if(CTS::BC == "FREE"){
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
        U1[0] = U1[1];
        U1[N+1] = U1[N];

        U2[0] = U2[1];
        U2[N+1] = U2[N];

        U3[0] = U3[1];
        U3[N+1] = U3[N];

    }
    else if(CTS::BC == "REFLECTIVE"){
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
        U1[0] = U1[1];
        U1[N+1] = U1[N];

        // Conserved variable-2
        U2[0] = -U2[1];
        U2[N+1] = -U2[N];

        // Conserved variable-3
        U3[0] = U3[1];
        U3[N+1] = U3[N];

    }
    else{
        cout << "---LOG--- Please enter correct boundary conditions---" << endl;
    }
}



