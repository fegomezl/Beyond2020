#include "header.h"

//Dimensions of the mesh
double Rmin;
double Zmin;
double Rmax;
double Zmax;

//State of the simulation
double mid;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);
    int refinements_init = 0;
    int refinements_last = 0;
    int DeltaT_init = 0;
    int DeltaT_jumpsize = 0;
    int DeltaT_jumps = 0;
    int state = 0;

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
    args.AddOption(&state, "-state", "--state",
                   "State of the simulation: \n"
                   "1.Onephase stable point.\n" 
                   "2.Twophase stable point.\n"
                   "3.Mid stable point.");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&refinements_init, "-r_i", "--refinements_init",
                  "Number of total initial uniform refinements.");
    args.AddOption(&refinements_last, "-r_l", "--refinements_last",
                  "Number of total last uniform refinements.");
    args.AddOption(&DeltaT_init, "-DT_i", "--DeltaT_init",
                   "Initial temperature interface interval (10^(-n)).");
    args.AddOption(&DeltaT_jumpsize, "-DT_js", "--DeltaT_jumpsize",
                   "Jump size of the temperature interface interval (10^(-n)).");
    args.AddOption(&DeltaT_jumps, "-DT_j", "--DeltaT_jumps",
                   "Jumps of the temperature interface interval (10^(-n)).");
    args.AddOption(&config.T_f, "-T_f", "--temperature_fusion",
                   "Fusion temperature of the material.");
    args.AddOption(&config.c_s, "-c_s", "--c_s",
                   "Solid volumetric heat capacity.");
    args.AddOption(&config.c_l, "-c_l", "--c_l",
                   "Liquid volumetric heat capacity.");
    args.AddOption(&config.k_s, "-k_s", "--k_s",
                   "Solid thermal conductivity.");
    args.AddOption(&config.k_l, "-k_l", "--k_l",
                   "Liquid thermal conductivity.");
    args.AddOption(&config.L, "-L", "--L",
                   "Volumetric latent heat.");
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

    if (state == 1)
        mid = Zmax/2;
    else if (state == 2)
        mid = Zmax*config.k_s/(config.k_s + config.k_l);
    else 
        mid = Zmax*(config.k_s/(config.k_s + config.k_l) + 0.5)/2.;

    for (config.refinements = refinements_init; config.refinements <= refinements_last; config.refinements++){
        for (config.nDeltaT = DeltaT_init; config.nDeltaT <= DeltaT_init + DeltaT_jumpsize*DeltaT_jumps; config.nDeltaT += DeltaT_jumpsize){
            //Run the program for different refinements
            config.invDeltaT = pow(10, config.nDeltaT);
            Artic_sea artic_sea(config);
            artic_sea.run(mesh_file);
        }
    }

    MPI_Finalize();

    return 0;
}
