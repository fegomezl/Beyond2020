#include "header.h"

//Dimensions of the mesh
double Rmin;
double Zmin;
double Rmax;
double Zmax;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);
    int nDeltaT = 0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&Rmin, "-Rmin", "--Rmin",
                   "Minimum R border.");
    args.AddOption(&Rmax, "-Rmax", "--Rmax",
                   "Maximum R border.");
    args.AddOption(&Zmin, "-Zmin", "--Zmin",
                   "Minimum Z border.");
    args.AddOption(&Zmax, "-Zmax", "--Zmax",
                   "Maximum Z border.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                   "Number of total uniform refinements.");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.ode_solver_type, "-ode", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.reltol, "-reltol", "--tolrelativaSUNDIALS",
                   "Tolerancia relativa de SUNDIALS solvers");
    args.AddOption(&config.abstol, "-abstol", "--tolabsolutaSUNDIALS",
                   "Tolerancia absoluta de S");
    args.AddOption(&config.T_f, "-T_f", "--temperature_fusion",
                   "Fusion temperature of the material.");
    args.AddOption(&nDeltaT, "-DT", "--DeltaT",
                   "Temperature interface interval (10^(-n)).");
    args.AddOption(&config.c_l, "-c_l", "--c_l",
                   "Liquid volumetric heat capacity.");
    args.AddOption(&config.c_s, "-c_s", "--c_s",
                   "Solid volumetric heat capacity.");
    args.AddOption(&config.k_l, "-k_l", "--k_l",
                   "Liquid thermal conductivity.");
    args.AddOption(&config.k_s, "-k_s", "--k_s",
                   "Solid thermal conductivity.");
    args.AddOption(&config.L, "-L", "--L",
                   "Volumetric latent heat.");

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "a\n";
    MPI_Barrier(MPI_COMM_WORLD);

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    //Run the program for different refinements
    config.invDeltaT = pow(10, nDeltaT);
    Artic_sea artic_sea(config);
    artic_sea.run(mesh_file);

    MPI_Finalize();

    return 0;
}
