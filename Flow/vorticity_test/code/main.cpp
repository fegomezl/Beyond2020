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
    int nDeltaT, nEpsilon_eta, nEpsilon_r;

    //Make -program parameters readeable in execution
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

    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&config.refinements, "-r", "--refinements",
                   "Number of total uniform refinements");

    args.AddOption(&config.T_f, "-T_f", "--temperature_fusion",
                   "Fusion temperature of the material.");
    args.AddOption(&nDeltaT, "-DT", "--DeltaT",
                   "Temperature interface interval (10^(-n)).");
    args.AddOption(&nEpsilon_eta, "-e_eta", "--epsilon_eta",
                   "Value of constatn epsilon for (1-phi)^2/(phi^3 + epsilon)(10^(-n)).");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    const char *petscrc_file = "settings/petsc_settings";
    MFEMInitializePetsc(NULL, NULL, petscrc_file, NULL);

    //Run the program
    {
        tic();
        config.invDeltaT = pow(10, nDeltaT);
        config.epsilon_eta = pow(10, -nEpsilon_eta);
        config.last = true;
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MFEMFinalizePetsc();

    MPI_Finalize();

    return 0;
}
