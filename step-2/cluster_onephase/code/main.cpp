#include "header.h"

double Rmin;
double Zmin;
double Rmax;
double Zmax;
double alpha;
int  Mterms=6;
int  Nterms=6;
std::vector<double> Coeficients(Mterms*Nterms ,0.0);

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&Rmin, "-Rmin", "--Rmin",
                   "Minimum R border");
    args.AddOption(&Rmax, "-Rmax", "--Rmax",
                   "Maximum R border");
    args.AddOption(&Zmin, "-Zmin", "--Zmin",
                   "Minimum Z boder");
    args.AddOption(&Zmax, "-Zmax", "--Zmax",
                   "Maximum Z boder");
    args.AddOption(&config.order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&config.refinements, "-r", "--refinements",
                  "Number of total uniform refinements.");
    args.AddOption(&alpha, "-a", "--alpha",
                   "Alpha coefficient for the material.");
    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&config.ode_solver_type, "-ode", "--ode_solver",
                   "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3, \n"
                   "            11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
    args.AddOption(&config.reltol, "-reltol", "--tolrelativaSUNDIALS",
                   "Tolerancia relativa de SUNDIALS solvers");
    args.AddOption(&config.abstol, "-abstol", "--tolabsolutaSUNDIALS",
                   "Tolerancia absoluta de S");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    //Calculate coeficients of the exact solution
    Calc_Coe(Zmax, Rmax, Coeficients);

    int refinements = config.refinements;
    int N = 2;
    vector<double> Data(3*(refinements + 1), 0.);

    //Run the program for different refinements
    for (config.refinements = 0; config.refinements <= refinements; config.refinements++){
        double h_size = 0, total_time = 0, mean_error = 0, null = 0;
        for (int ii = 0; ii <= N; ii++){
            if (ii != N){
                config.print = false;
                Artic_sea artic_sea(config);
                artic_sea.run(mesh_file, h_size, total_time, mean_error);
            } else {
                config.print = true;
                Artic_sea artic_sea(config);
                artic_sea.run(mesh_file, null, null, null);
            }
        }

        Data[3*config.refinements] = h_size/N;
        Data[3*config.refinements + 1] = total_time/N;
        Data[3*config.refinements + 2] = mean_error/N;
    }

    if (config.master){
        string fname = "results/global_" + to_string(config.ode_solver_type) + "_" + to_string(config.order) + ".txt";
        ofstream output;
        output.precision(4);
        output.open(fname, ios::trunc);
        for (int ii = 0; ii <= refinements; ii++)
            output << left << setw(12)
                   << Data[3*ii] << setw(12)
                   << Data[3*ii + 1] << setw(12)
                   << Data[3*ii + 2]  << "\n";
        output.close();
    }

    MPI_Finalize();

    return 0;
}
