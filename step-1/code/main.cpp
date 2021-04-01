#include "header.h"

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    bool master = (pid == 0);

    //Define program paramenters
    int order = 1;
    int refinements = 6;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&refinements, "-r", "--refinements",
                  "Number of total uniform refinements");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
      if (master) args.PrintUsage(cout);
      MPI_Finalize();
      return 1;
    }
    if (master) args.PrintOptions(cout);

    Artic_sea artic_sea(master, order, refinements);
    artic_sea.run();

    MPI_Finalize();

    return 0;
}
