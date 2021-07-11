#include "header.h"

//Dimensions of the mesh
double height;
double out_rad;
double int_rad;

//Size of the BC border
double border;
double InvR;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);
    int nDeltaT, nCold_porosity, nInvR;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&height, "-he", "--height",
                   "Height of the container.");
    args.AddOption(&out_rad, "-o-r", "--outer-radius",
                   "Outer radius of the container.");
    args.AddOption(&int_rad, "-i-r", "--internal-radius",
                   "Internal radius of the container.");

    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&config.refinements, "-r", "--refinements",
                   "Number of total uniform refinements");
    args.AddOption(&border, "-b", "--border",
                   "Border size for the BCs in fraction (1/x).");

    args.AddOption(&config.T_f, "-T_f", "--temperature_fusion",
                   "Fusion temperature of the material.");
    args.AddOption(&nDeltaT, "-DT", "--DeltaT",
                   "Temperature interface interval (10^(-n)).");
    args.AddOption(&config.viscosity, "-v", "--viscosity",
                   "Kinematic viscosity of the material.");
    args.AddOption(&nCold_porosity, "-ct", "--cold_porosity",
                   "Value of the porosity on the solid domain (10^(n)).");
    args.AddOption(&nInvR, "-ir", "--inv_r",
                   "Value of constant m for 1/(r + m) (10^(n)).");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    //Run the program for different refinements
    int total_refinements = config.refinements;
    for (int ii = total_refinements; ii <= total_refinements; ii++){
        config.invDeltaT = pow(10, nDeltaT);
        config.cold_porosity = pow(10, -nCold_porosity);
        InvR = pow(10, -nInvR);
        config.last = ((config.refinements = ii) == total_refinements);
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
