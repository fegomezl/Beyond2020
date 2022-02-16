#include "header.h"

//Dimensions of the mesh
double R;
double Z;

//Size of the indetermination window in jump functions
double JumpScale;
double Epsilon;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);
    int nJumpScale = 0;
    double T_l, T_s;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&R, "-R", "--Radius",
                   "Radius of the cilinder.");
    args.AddOption(&Z, "-Z", "--Height",
                   "Height of the cilinder.");

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
    args.AddOption(&config.abstol_conduction, "-abstol_c", "--tolabsoluteConduction",
                   "Absolute tolerance of Conduction.");
    args.AddOption(&config.reltol_conduction, "-reltol_c", "--tolrelativeConduction",
                   "Relative tolerance of Conduction.");
    args.AddOption(&config.iter_conduction, "-iter_c", "--iterationsConduction",
                   "Iterations of Conduction.");
    args.AddOption(&config.abstol_sundials, "-abstol_s", "--tolabsoluteSUNDIALS",
                   "Absolute tolerance of SUNDIALS.");
    args.AddOption(&config.reltol_sundials, "-reltol_s", "--tolrelativeSUNDIALS",
                   "Relative tolerance of SUNDIALS.");

    args.AddOption(&nJumpScale, "-nJ", "--nJumpScale",
                   "Inverse of the lenght of state indetermination in jump functions (10^(-n)).");
    args.AddOption(&config.a_l, "-a_l", "--a_l",
                   "Liquid heat diffusion constant.");
    args.AddOption(&config.a_s, "-a_s", "--a_s",
                   "Solid heat diffusion constant.");

    args.AddOption(&T_l, "-T_l", "--T_liquid",
                   "Liquid temperature interval.");
    args.AddOption(&T_s, "-T_s", "--T_solid",
                   "Solid temperature interval.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    //Run the program
    {
        tic();
        JumpScale = pow(10, nJumpScale);
        Epsilon = pow(10, -nJumpScale);
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
