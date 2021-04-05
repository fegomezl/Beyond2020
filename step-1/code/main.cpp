#include "header.h"

double height = 20;
double int_rad = 5;
double out_rad = 20;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    bool master = (pid == 0);

    //Define program paramenters
    const char *mesh_file = "data/transfinite_cylinder.msh";
    int order = 1;
    int refinements = 0;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&refinements, "-r", "--refinements",
                  "Number of total uniform refinements");
    args.AddOption(&height, "-h", "--height",
                   "Height of the container.");
    args.AddOption(&int_rad, "-i-r", "--internal-radius",
                   "Internal radius of the container.");
    args.AddOption(&out_rad, "-o-r", "--outer-radius",
                   "Outer radius of the container.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (master) args.PrintOptions(cout);

    //Run the program for different refinements
    for (int ii = 0; ii <= refinements; ii++){
        Artic_sea artic_sea(master, order, ii);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
