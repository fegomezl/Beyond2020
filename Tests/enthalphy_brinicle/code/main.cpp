#include "header.h"

double R;
double Z;
double R_in;
double Z_out;

double LenghtScale;  
double TimeScale;    
double JumpScale;
double Epsilon;

double InflowVelocity;
double InflowFlux;
double InitialTemperature;       
double InflowTemperature;        
double InitialSalinity;          
double InflowSalinity;           
double NucleationLength;         
double NucleationHeight;         
double NucleationTemperature;    
double NucleationSalinity;       

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config((pid == 0), nproc);
    int nEpsilon = 0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&R, "-R", "--Radius",
                   "Radius of the cilinder.");
    args.AddOption(&Z, "-Z", "--Height",
                   "Height of the cilinder.");
    args.AddOption(&R_in, "-Ri", "--Radius_in",
                   "Radius of inflow.");
    args.AddOption(&Z_out, "-Zo", "--Height_out",
                   "Height of outflow.");

    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    args.AddOption(&LenghtScale, "-LS", "--lenght_scale",
                   "Scaling of the lenght dimension (currently in mm).");
    args.AddOption(&TimeScale, "-TS", "--time_scale",
                   "Scaling of the time dimension (currently in min).");

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

    args.AddOption(&nEpsilon, "-nE", "--nEpsilon",
                   "Lenght of state indetermination in heaviside functions (10^(-n)).");

    args.AddOption(&InflowVelocity, "-V", "--velocity",
                   "Velocity of the inflow.");
    args.AddOption(&InitialTemperature, "-To", "--initial_temperature",
                   "Initial temperature of the domain.");
    args.AddOption(&InflowTemperature, "-Ti", "--inflow_temperature",
                   "Temperature of the inflow.");
    args.AddOption(&InitialSalinity, "-So", "--initial_salinity",
                   "Initial salinity of the domain.");
    args.AddOption(&InflowSalinity, "-Si", "--inflow_salinity",
                   "Salinity of the inflow.");
    args.AddOption(&NucleationLength, "-Nl", "--nucleation_lenght",
                   "Lenght of the nucleation point.");
    args.AddOption(&NucleationHeight, "-Nh", "--nucleation_height",
                   "Height of the nucleation point.");
    args.AddOption(&NucleationTemperature, "-Nt", "--nucleation_temperature",
                   "Temperature of the nucleation point.");
    args.AddOption(&NucleationSalinity, "-Ns", "--nucleation_salinity",
                   "Salinity of the nucleation point.");

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
        InflowFlux = 0.25*InflowVelocity*pow(R_in, -2);
        Epsilon = pow(10, -nEpsilon);
        JumpScale = pow(10, nEpsilon);
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

return 0;
}
