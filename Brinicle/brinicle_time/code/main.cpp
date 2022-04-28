#include "header.h"

/****
 * All the nomeclarure is explained in the header
 ****/

double RMin, RMax, ZMin, ZMax;
double RIn, ZOut;
double Epsilon, EpsilonInv;

double InflowVelocity;   
double InitialTemperature;
double InflowTemperature;
double InitialSalinity;
double InflowSalinity;
double NucleationLength;
double NucleationHeight;
double NucleationTemperature;
double NucleationSalinity;

double InflowFlux;                    
double RelaxationTime = 0.08;
double RenormalizationScale = 1.E+9;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config(pid, nproc);
    int rescale = 0;
    int nEpsilon = 0;
    int restart = 0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&RMin, "-Rmin", "--Rmin",
                   "Minimum R border.");
    args.AddOption(&RMax, "-Rmax", "--Rmax",
                   "Maximum R border.");
    args.AddOption(&ZMin, "-Zmin", "--Zmin",
                   "Minimum Z border.");
    args.AddOption(&ZMax, "-Zmax", "--Zmax",
                   "Maximum Z border.");
    args.AddOption(&RIn, "-Li", "--L_in",
                   "Inflow window size.");
    args.AddOption(&ZOut, "-Lo", "--L_out",
                   "Outflow window size.");

    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&rescale, "-rc", "--rescale",
                   "If the simulation rescales the stream (1) or not (0).");

    args.AddOption(&config.refinements, "-ref", "--refinements",
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
    args.AddOption(&nEpsilon, "-eps", "--epsilon",
                   "Epsilon constant for heaviside functions (10^(-n)).");

    args.AddOption(&InflowVelocity, "-v", "--vel",
                   "Inflow velocity.");
    args.AddOption(&InitialTemperature, "-Ti", "--Theta_in",
                   "Initial temperature.");
    args.AddOption(&InflowTemperature, "-To", "--Theta_out",
                   "Inflow temperature.");
    args.AddOption(&InitialSalinity, "-Si", "--Phi_in",
                   "Initial salinity.");
    args.AddOption(&InflowSalinity, "-So", "--Phi_out",
                   "Inflow salinity.");
    args.AddOption(&NucleationLength, "-nl", "--n_length",
                   "Nucleation length.");
    args.AddOption(&NucleationHeight, "-nh", "--n_heigth",
                   "Nucleation height.");
    args.AddOption(&NucleationTemperature, "-Tn", "--Theta_n",
                   "Nucleation temperature.");
    args.AddOption(&NucleationSalinity, "-Sn", "--Phi_n",
                   "Nucleation salinity.");

    args.AddOption(&restart, "-r", "--restart",
                   "If the simulation restarts (1) or not (0).");
    args.AddOption(&config.t_init, "-t_i", "--t_init",
                   "Start time of restart.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    {
        InflowFlux = 0.25*InflowVelocity*pow(RIn, 2);

        config.rescale = (rescale == 1);

        Epsilon = pow(10, -nEpsilon); 
        EpsilonInv = pow(10, nEpsilon); 

        config.restart = (restart == 1);
        config.t_init = config.restart ? config.t_init : 0.;
        config.t_final += config.t_init;

        tic();
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
