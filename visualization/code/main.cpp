#include "header.h"

double T_f;       //Fusion temperature
double T_i;      //Initial temperature

double alpha_l;  //Liquid thermal conduction
double alpha_s; //Solid thermal conduction
double lambda; //s(t) = sqrt(4*lamda*(alpha_s+alpha_l)*t)

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&T_f, "-T_f", "--temperature_fusion",
                   "Fusion Temperature of the material.");
    args.AddOption(&T_i, "-T_i", "--temperature_initial",
                   "Initial temperature of the material.");
    args.AddOption(&alpha_l, "-a_l", "--alpha_liquid",
                   "Alpha coefficient for liquid phase.");
    args.AddOption(&alpha_s, "-a_s", "--alpha_solid",
                   "Alpha coefficient for solid phase.");
    args.AddOption(&lambda, "-l", "--lambda",
                   "Lambda constant.");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                  "Number of total uniform refinements.");
    args.AddOption(&config.dt, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_init, "-t_i", "--t_init",
                   "Initial time.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");

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
