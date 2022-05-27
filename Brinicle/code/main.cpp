#include "header.h"

/*
 * All the nomeclarure is explained in the header
 *
 */

double L_ref;
double V_ref;
double t_ref;
double T0_ref;
double T_ref;
double S0_ref;
double S_ref;

double R, Z;
double RInflow;
double Epsilon, EpsilonInv;

double FluxRate;
double InflowFlux;                    
double NucleationLength;
double NucleationHeight;
double InitialTemperature;
double InflowTemperature;
double NucleationTemperature;
double InitialSalinity;
double InflowSalinity;
double NucleationSalinity;

int main(int argc, char *argv[]){
    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    //Define program paramenters
    const char *mesh_file;
    Config config(pid, nproc);
    int nEpsilon = 0;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&R, "-R", "--R",
                   "Radius of the domain cylinder.");
    args.AddOption(&Z, "-Z", "--Z",
                   "Lenght of the domain cylinder.");
    args.AddOption(&RInflow, "-R_in", "--R_inflow",
                   "Inflow window radius.");

    args.AddOption(&config.dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&config.t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&config.vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");

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

    args.AddOption(&FluxRate, "-q", "--Q",
                   "Flux of the inflow boundary.");
    args.AddOption(&NucleationLength, "-nl", "--n_length",
                   "Nucleation length.");
    args.AddOption(&NucleationHeight, "-nh", "--n_heigth",
                   "Nucleation height.");
    args.AddOption(&InitialTemperature, "-Ti", "--Theta_in",
                   "Initial temperature.");
    args.AddOption(&InflowTemperature, "-To", "--Theta_out",
                   "Inflow temperature.");
    args.AddOption(&NucleationTemperature, "-Tn", "--Theta_n",
                   "Nucleation temperature.");
    args.AddOption(&InitialSalinity, "-Si", "--Phi_in",
                   "Initial salinity.");
    args.AddOption(&InflowSalinity, "-So", "--Phi_out",
                   "Inflow salinity.");
    args.AddOption(&NucleationSalinity, "-Sn", "--Phi_n",
                   "Nucleation salinity.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
        if (config.master) args.PrintUsage(cout);
        MPI_Finalize();
        return 1;
    }
    if (config.master) args.PrintOptions(cout);

    {
        L_ref = RInflow;
        V_ref = FluxRate*pow(L_ref, -2);
        t_ref = L_ref/V_ref;
        T0_ref = InflowTemperature;
        T_ref = InitialTemperature-InflowTemperature;
        S0_ref = InitialSalinity;
        S_ref = InflowSalinity-InitialSalinity;

        config.dt_init /= t_ref;
        config.t_final /= t_ref; 

        R /= L_ref;
        Z /= L_ref;
        RInflow /= L_ref;

        FluxRate = 1.;
        InflowFlux = 0.5*M_1_PI;
        NucleationLength /= L_ref;
        NucleationHeight /= L_ref;

        InitialTemperature = (InitialTemperature-T0_ref)/T_ref;
        InflowTemperature = (InflowTemperature-T0_ref)/T_ref;
        NucleationTemperature = (NucleationTemperature-T0_ref)/T_ref;

        InitialSalinity = (InitialSalinity-S0_ref)/S_ref;
        InflowSalinity = (InflowSalinity-S0_ref)/S_ref;
        NucleationSalinity = (NucleationSalinity-S0_ref)/S_ref;

        Epsilon = pow(10, -nEpsilon); 
        EpsilonInv = pow(10, nEpsilon); 

        if (config.master)
            cout << "\nAdimentional numbers:\n"
                 << "Reynolds number: " << L_ref*V_ref/constants.Viscosity << "\n" 
                 << "Froude number: " << V_ref*pow(constants.Gravity*L_ref, -0.5) << "\n" 
                 << "Peclet number(T_l): " << L_ref*V_ref/constants.TemperatureDiffusion_l << "\n" 
                 << "Peclet number(T_s): " << L_ref*V_ref/constants.TemperatureDiffusion_s << "\n" 
                 << "Peclet number(S_l): " << L_ref*V_ref/constants.SalinityDiffusion_l << "\n" 
                 << "Peclet number(S_s): " << L_ref*V_ref/constants.SalinityDiffusion_s << "\n" 
                 << "Stefan number: " << T_ref/constants.Stefan << "\n\n";

        tic();
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
