#include "header.h"

double Rmin = 2;
double Zmin = 0;
double Rmax = 10;
double Zmax = 10;
double T_f = 0;
double alpha_l = 7.8;
double alpha_s = 70.8;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);


    config.order = 2;
    config.refinements = 3;
    config.dt_init = 0.001;
    config.t_final = 0.5;
    config.vis_steps = 10;
    config.ode_solver_type = 12;
    config.reltol = config.abstol = 0.00001;

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
    args.AddOption(&alpha_l, "-a_l", "--alpha_liquid",
                   "Alpha coefficient for liquid phase.");
    args.AddOption(&alpha_s, "-a_s", "--alpha_solid",
                   "Alpha coefficient for solid phase.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&config.ode_solver_type, "-ode", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.reltol, "-restol", "--tolrelativaSUNDIALS",
                   "Initial time step.");
    args.AddOption(&config.abstol, "abstol", "--tolabsolutaSUNDIALS",
                   "Initial time step.");

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
