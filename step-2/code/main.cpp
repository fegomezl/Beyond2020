#include "header.h"

double T_f = 0.;       //Fusion temperature
double T_i = 10.;      //Initial temperature

double alpha_s = 70.8; //Solid thermal conduction
double alpha_l = 7.8;  //Liquid thermal conduction
double lamda = 0.0721; //s(t) = sqrt(4*lamda*(alpha_s+alpha_l)*t)

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config(pid == 0);
    config.serial_refinements = 0;
    config.ode_solver_type = 3;
    config.dt_init = 0.001;
    config.t_final = 0.02;
    config.vis_steps = 10;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.serial_refinements, "-s_r", "--serial_refinements",
                  "Number of serial uniform refinements.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                  "Number of total uniform refinements.");
    args.AddOption(&config.ode_solver_type, "-s", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&T_f, "-T_f", "--Temperature_fusion",
                   "Fusion Temperature of the material.");
    args.AddOption(&T_i, "-T_i", "--Temperature_initial",
                   "Initial temperature of the material.");
    args.AddOption(&alpha_l, "-a_l", "--alpha_liquid",
                   "Alpha coefficient for liquid phase.");
    args.AddOption(&alpha_s, "-a_s", "--alpha_solid",
                   "Alpha coefficient for solid phase.");
    args.AddOption(&lamda, "-l", "--lamda",
                   "Lamda constant.");
    args.AddOption(&config.vis_steps, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");

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
