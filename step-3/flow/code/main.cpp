#include "header.h"

//Dimensions of the mesh
double height;
double out_rad;
double int_rad;
double k1 = 1;
double k2 = 1;

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
    for (int ii = 0; ii <= total_refinements; ii++){
        config.last = ((config.refinements = ii) == total_refinements);
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
