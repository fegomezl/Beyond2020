#include "header.h"

double Rmin;
double Zmin;
double Rmax;
double Zmax;
double T_f;
double c_s = 1.878, c_l = 4.219;
double k_s = 133.2, k_l = 33.6;
double L = 302.3;

double DeltaT = 0.0001;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&Rmin, "-Rmin", "--Rmin",
                   "Minimum R border");
    args.AddOption(&Rmax, "-Rmax", "--Rmax",
                   "Maximum R border");
    args.AddOption(&Zmin, "-Zmin", "--Zmin",
                   "Minimum Z border");
    args.AddOption(&Zmax, "-Zmax", "--Zmax",
                   "Maximum Z boder");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                  "Number of total uniform refinements.");
    args.AddOption(&T_f, "-T_f", "--temperature_fusion",
                   "Fusion Temperature of the material.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&config.ode_solver_type, "-ode", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.reltol, "-reltol", "--tolrelativaSUNDIALS",
                   "Tolerancia relativa de SUNDIALS solvers");
    args.AddOption(&config.abstol, "-abstol", "--tolabsolutaSUNDIALS",
                   "Tolerancia absoluta de S");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    //Run the program for different refinements
    Artic_sea artic_sea(config);
    artic_sea.run(mesh_file);

    MPI_Finalize();

    return 0;
}
