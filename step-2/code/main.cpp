#include "header.h"

double out_rad;
double int_rad;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    bool master = (pid == 0);

    //Define program paramenters
    const char *mesh_file;
    Config config;
    config.master = (pid == 0);
    config.serial_refinements = 0;
    config.ode_solver_type = 3;
    config.t_init = 0.;
    config.t_final = 1.;
    config.dt_init = 1.0e-2;
    config.alpha = 1.0e-2;
    config.kappa = 0.5;
    config.vis_steps = 10;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&out_rad, "-o_r", "--outer_radius",
                   "Outer radius of the container.");
    args.AddOption(&int_rad, "-i_r", "--internal_radius",
                   "Internal radius of the container.");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.serial_refinements, "-s_r", "--serial_refinements",
                  "Number of serial uniform refinements.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                  "Number of total uniform refinements.");
    args.AddOption(&config.ode_solver_type, "-s", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.t_init, "-t_i", "--t_initial",
                   "Starting time.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.alpha, "-a", "--alpha",
                   "Alpha coefficient.");
    args.AddOption(&config.kappa, "-k", "--kappa",
                   "Kappa coefficient offset.");
    args.AddOption(&config.vis_steps, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (master) args.PrintOptions(cout);

    //Run the program for different refinements
    Artic_sea artic_sea(config);
    artic_sea.run(mesh_file);

    MPI_Finalize();

    return 0;
}
