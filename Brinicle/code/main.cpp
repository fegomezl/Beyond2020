#include "header.h"

Constants constants;

double L_ref;
double t_ref;
double T0_ref;
double T_ref;
double S0_ref;
double S_ref;

double R, Z;
double RInflow;
double Epsilon, EpsilonInv;

double FluxRate;
double NucleationLength;
double NucleationHeight;
double InitialTemperature;
double InflowTemperature;
double NucleationTemperature;
double ZeroTemperature;
double InitialSalinity;
double InflowSalinity;
double NucleationSalinity;
double ZeroSalinity;

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
        t_ref = 2*M_PI*pow(RInflow, 3)/FluxRate;
        T0_ref = InitialTemperature;
        T_ref = InflowTemperature-InitialTemperature;
        S0_ref = InitialSalinity;
        S_ref = InflowSalinity-InitialSalinity;

        config.dt_init /= t_ref;
        config.t_final /= t_ref; 

        R /= L_ref;
        Z /= L_ref;
        RInflow = 1.;

        NucleationLength /= L_ref;
        NucleationHeight /= L_ref;

        InitialTemperature = 0.;
        InflowTemperature = 1.;
        NucleationTemperature = (NucleationTemperature-T0_ref)/T_ref;
        ZeroTemperature = T0_ref/T_ref;

        InitialSalinity = 0.;
        InflowSalinity = 1.;
        NucleationSalinity = (NucleationSalinity-S0_ref)/S_ref;
        ZeroSalinity = S0_ref/S_ref;

        constants.FusionPoint_a *= S_ref/T_ref;
        constants.FusionPoint_b *= pow(S_ref, 3)/T_ref;

        constants.Stefan /= abs(T_ref);

        constants.TemperatureDiffusion_l *= t_ref*pow(L_ref, -2);
        constants.TemperatureDiffusion_s *= t_ref*pow(L_ref, -2);
        constants.SalinityDiffusion_l *= t_ref*pow(L_ref, -2);
        constants.SalinityDiffusion_s *= t_ref*pow(L_ref, -2);

        constants.Density_a0 *= S_ref;
        constants.Density_a1 *= S_ref*T_ref;
        constants.Density_a2 *= S_ref*pow(T_ref, 2);
        constants.Density_a3 *= S_ref*pow(T_ref, 3);
        constants.Density_a4 *= S_ref*pow(T_ref, 4);
        constants.Density_b0 *= pow(abs(S_ref), 1.5);
        constants.Density_b1 *= pow(abs(S_ref), 1.5)*T_ref;
        constants.Density_b2 *= pow(abs(S_ref), 1.5)*pow(T_ref, 2);
        constants.Density_c0 *= pow(S_ref, 2);

        constants.BuoyancyCoefficient *= t_ref*L_ref;

        Epsilon = pow(10, -nEpsilon); 
        EpsilonInv = pow(10, nEpsilon); 

        if (config.master){
            cout << "\na:\n"
                 << "T: " << T0_ref << "\n"
                 << "S: " << S0_ref << "\n";
            cout << "Fusion Point:\n"
                 << "a: " << constants.FusionPoint_a << "\n"
                 << "b: " << constants.FusionPoint_b << "\n"
                 << "Density:\n"
                 << "a0: " << constants.Density_a0 << "\n"
                 << "a1: " << constants.Density_a1 << "\n"
                 << "a2: " << constants.Density_a2 << "\n"
                 << "a3: " << constants.Density_a3 << "\n"
                 << "a4: " << constants.Density_a4 << "\n"
                 << "b0: " << constants.Density_b0 << "\n"
                 << "b1: " << constants.Density_b1 << "\n"
                 << "b2: " << constants.Density_b2 << "\n"
                 << "c0: " << constants.Density_c0 << "\n";
            cout << "Adimentional numbers:\n"
                 << "Reynolds number: "    << pow(RInflow, 2)/(t_ref*constants.Viscosity) << "\n" 
                 << "Froude number: "      << pow(RInflow/constants.Gravity, 0.5)/t_ref << "\n" 
                 << "Peclet number(T_l): " << 1/constants.TemperatureDiffusion_l << "\n" 
                 << "Peclet number(T_s): " << 1/constants.TemperatureDiffusion_s << "\n" 
                 << "Peclet number(S_l): " << 1/constants.SalinityDiffusion_l << "\n" 
                 << "Peclet number(S_s): " << 1/constants.SalinityDiffusion_s << "\n" 
                 << "Stefan number: " << 1/constants.Stefan << "\n";
        }

        tic();
        Artic_sea artic_sea(config);
        artic_sea.run(mesh_file);
    }

    MPI_Finalize();

    return 0;
}
