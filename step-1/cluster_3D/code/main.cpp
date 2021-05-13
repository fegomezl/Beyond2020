#include "header.h"

double height;
double out_rad;
double int_rad;

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

    int refinements = config.refinements;
    int N = 2;
    vector<double> Data(3*(refinements + 1), 0.);

    //Run the program for different refinements
    for (config.refinements = 0; config.refinements <= refinements; config.refinements++){
        double h_size = 0, total_time = 0, mean_error = 0;
        for (int ii = 0; ii <= N; ii++){
            Artic_sea artic_sea(config);
            artic_sea.run(mesh_file, h_size, total_time, mean_error);
        }

        Data[3*config.refinements] = h_size/N;
        Data[3*config.refinements + 1] = total_time/N;
        Data[3*config.refinements + 2] = mean_error/N;
    }

    if (config.master){
        string fname = "results/global_" + to_string(config.order) + ".txt";
        ofstream output;
        output.precision(4);
        output.open(fname, ios::trunc);
        for (int ii = 0; ii <= refinements; ii++)
            output << left << setw(12)
                   << Data[3*ii] << setw(12)
                   << Data[3*ii + 1] << setw(12)
                   << Data[3*ii + 2] << "\n";
        output.close();
    }

    MPI_Finalize();

    return 0;
}
