#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

//Main constants for the simulation
struct Constants{ 
    /****
     * Physical properties of the system
     *
     * Density of liquid phase (rho_l): 1028.3 kg/m^3
     * Specific heat capacity of liquid phase (c_l): 3.77 kJ/(kg*°C)
     * Thermal conductivity of liquid phase (k_l): 0.5433 W/(m*°C)
     * Salt diffusivity of liquid phase (d_l): 0.001 mm^2/s
     *
     * Density of solid phase (rho_s): 934.4 kg/m^3
     * Specific heat capacity of solid phase (c_s): 2.12 kJ/(kg*°C)
     * Thermal conductivity of solid phase (k_s): 2.2467 W/(m*°C)
     * Salt diffusivity of solid phase (d_s): 0 mm^2/s
     *
     * Kinematic viscosity (nu): 6.8 mm^2/s
     * Latent heat (L): 322.7 kJ/kg
     * Gravity (g): 9.8 m/s^2
     ****/

    /****
     * Coefficients for the equation of the fusion 
     * temperature T_f in terms of salinity S:
     * T_f = -(a*S + b *S^3)
     * in °C
     ****/
    double FusionPoint_a = 6.04E-1;
    double FusionPoint_b = 5.81E-4;

    /****
     * Coefficients for the mass term of the temperature
     * equation given by specific heat capacity over 
     * latent heat (c/L) in 1/°C
     ****/ 
    double TemperatureMass_l = 1.17E-2;      //Liquid;
    double TemperatureMass_s = 6.57E-3;      //Solid;

    /****
     * Coefficients for the diffusion term of the temperature 
     * equation given by thermal conductivity over volumetric 
     * latent heat (k/(rho*L)) in  mm^2/min
     ****/ 
    double TemperatureDiffusion_l = 9.82E-2;      //Liquid;
    double TemperatureDiffusion_s = 4.47E-1;      //Solid;

    /****
     * Coefficients for the diffusion term of the salinity
     * equation given by salt diffusivity (d) in mm^2/min
     ****/ 
    double SalinityDiffusion_l = 5.00E+0;       //Liquid
    double SalinityDiffusion_s = 1.00E-7;       //Solid

    /****
     * Coefficients for the relative density equation in terms of 
     * temperature T and salinity S given by:
     * Delta rho / rho = (a0 + a1*T + a2*T^3 + a3*T^3 + a4*T^4)*S
     *                 + (b0 + b1*T + b2*T^2)*S^1.5
     *                 + (c0)*S^2
     * which is dimentionless
     ****/
    double Density_a0 = 8.25E-3;
    double Density_a1 = -4.11E-5;
    double Density_a2 = 7.73E-7;
    double Density_a3 = -8.32E-9;
    double Density_a4 = 5.52E-11;
    double Density_b0 = -1.83E-4;
    double Density_b1 = 3.33E-6;
    double Density_b2 = -5.60E-8;
    double Density_c0 = 4.89E-5;

    /****
     * Coefficient for the buoyancy term(k) 
     * given by g/nu in 1/(mm*min)
     ****/ 
    double BuoyancyCoefficient = 8.647E+4;
};

//Main variables for the program
struct Config{
    //Initialization of the program variables 
    Config(int pid, int nproc);     

    //MPI variables
    bool master;
    int nproc;
    int pid;

    //Time and visualization variables
    double dt_init;
    double t_final;
    int vis_steps_max;
    bool rescale;

    //FEM variables 
    int refinements;
    int order;
    double reltol_conduction;
    double abstol_conduction;
    int iter_conduction;
    double reltol_sundials;
    double abstol_sundials;

    //Re-Initialization variables
    bool restart;
    double t_init;
};


//Solver for the temperature and salinity field
class Transport_Operator : public TimeDependentOperator{
    public:
        //Initialization of the solver
        Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X);

        //Update of the solver on each iteration
        void SetParameters(const BlockVector &X, const Vector &rVelocity);

        //Time-evolving functions
        virtual void Mult(const Vector &X, Vector &dX_dt) const;
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);

        virtual ~Transport_Operator();
    protected:
        //All 0-variables are related to temperature
        //All 1-variables are related to salinity

        //Global parameters
        Config config;

        //FEM parameters
        ParFiniteElementSpace &fespace_H1;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;
        Array<int> ess_bdr_0, ess_bdr_1;

        //Auxiliar grid functions
        ParGridFunction temperature, salinity, phase;
        ParGridFunction rvelocity;
        ParGridFunction heat_inertia;
        ParGridFunction heat_diffusivity;
        ParGridFunction salt_diffusivity;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient coeff_zero;

        ProductCoefficient coeff_rM;
        ProductCoefficient coeff_rD0; 
        ProductCoefficient coeff_rD1;

        VectorGridFunctionCoefficient coeff_rV;
        ScalarVectorProductCoefficient coeff_rMV;

        InnerProductCoefficient coeff_dPdT;
        InnerProductCoefficient coeff_dT_2;

        //System objects
        HypreParMatrix *M0, *M1;
        HypreParMatrix *M0_e, *M1_e; 
        HypreParMatrix *M0_o, *M1_o;
        HypreParMatrix *K0, *K1;
        HypreParMatrix *T0, *T1;
        HypreParMatrix *T0_e, *T1_e;

        HypreParVector *B0, *B1;
        HypreParVector B0_dt, B1_dt;
        mutable HypreParVector Z0, Z1;

        //Solver objects
        HyprePCG M0_solver, M1_solver;
        HyprePCG T0_solver, T1_solver;
        HypreBoomerAMG M0_prec, M1_prec;
        HypreBoomerAMG T0_prec, T1_prec;
};

//Solver for the velocity field
class Flow_Operator{
    public:
        //Initialization of the solver
        Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1);

        //Update of the solver on each iteration
        void SetParameters(const BlockVector &X);

        //Solution of the current system
        void Solve(BlockVector &Y, Vector &Velocity, Vector &rVelocity);

        ~Flow_Operator();
    protected:
        //All 0-variables are related to vorticity
        //All 1-variables are related to stream

        //Global parameters
        Config config;

        //FEM objects
        ParFiniteElementSpace &fespace_H1;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;
        Array<int> ess_bdr_0, ess_bdr_1;
        Array<int> ess_bdr_in, ess_bdr_out;
        Array<int> ess_bdr_closed_down, ess_bdr_closed_up;

        //Auxiliar grid functions
        ParGridFunction vorticity_boundary;
        ParGridFunction stream_boundary;
        ParGridFunction velocity;
        ParGridFunction rvelocity;
        ParGridFunction stream;
        ParGridFunction stream_gradient;

        ParGridFunction temperature;
        ParGridFunction salinity;
        ParGridFunction density;
        ParGridFunction density_dr;
        ParGridFunction impermeability;

        //System objects
        BlockVector B;
        HypreParVector *B0, *B1;

        HypreParMatrix *A00;
        HypreParMatrix *A01;
        HypreParMatrix *A10;
        HypreParMatrix *A11;
      
        //Coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient coeff_r_inv;
        VectorFunctionCoefficient coeff_r_inv_hat;
        MatrixFunctionCoefficient coeff_rot;

        FunctionCoefficient coeff_vorticity;
        FunctionCoefficient coeff_stream;
        FunctionCoefficient coeff_stream_in;
        FunctionCoefficient coeff_stream_out;
        ConstantCoefficient coeff_stream_closed_down;
        ConstantCoefficient coeff_stream_closed_up;

        //Interpolation objects
        ParDiscreteLinearOperator gradient;
};

//Main class of the program
class Artic_sea{
    public:
        //Initialization of the program
        Artic_sea(Config config);

        //Run the program
        void run(const char *mesh_file);

        ~Artic_sea();
    private:
        //Create the mesh and the FES
        void make_grid(const char *mesh_file);

        //Initialize the solvers and the variables
        void assemble_system();

        //Evolve the simulation one time step 
        void time_step();

        //Print the final results
        void output_results();

        //Global parameters
        Config config;

        //Simulation parameters
        int iteration;
        double t;
        double dt;
        bool last;
        int vis_iteration;
        int vis_steps;
        int vis_print;
        double total_time;

        //FEM objects
        ParMesh *pmesh;

        FiniteElementCollection *fec_H1;
        FiniteElementCollection *fec_ND;

        ParFiniteElementSpace *fespace_H1;
        ParFiniteElementSpace *fespace_ND;

        Array<int> block_offsets_H1;
        
        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size_H1;
        HYPRE_Int size_ND;

        //System objects
        ParGridFunction *temperature;
        ParGridFunction *salinity;
        ParGridFunction *phase;
        ParGridFunction *vorticity;
        ParGridFunction *stream;
        ParGridFunction *velocity;
        ParGridFunction *rvelocity;

        BlockVector X;
        BlockVector Y;
        HypreParVector *Velocity;
        HypreParVector *rVelocity;

        //Solvers
        Transport_Operator *transport_oper;
        Flow_Operator *flow_oper;

        //Time evolving operators
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Output gate
        ParaViewDataCollection *paraview_out;
};

//Constants associated with physical properties
const static Constants constants;

//Simulation parameters
extern double RMin, RMax, ZMin, ZMax;       //Size of the domain
extern double RIn, ZOut;                    //Size of the inflow and outflow
extern double Epsilon, EpsilonInv;          //Size of the indetermination window in heaviside functions 
                                            //(10̣^(-n) and 10^(n) respectively)

//Brinicle conditions
extern double InflowVelocity;           //Velocity of the inflow
extern double InitialTemperature;       //Initial temperature of the domain
extern double InflowTemperature;        //Temperature of the inflow
extern double InitialSalinity;          //Initial salinity of the domain
extern double InflowSalinity;           //Salinity of the inflow
extern double NucleationLength;         //Lenght of the nucleation point
extern double NucleationHeight;         //Height of the nucleation point
extern double NucleationTemperature;    //Temperature of the nucleation point
extern double NucleationSalinity;       //Salinity of the nucleation point
                                            
extern double InflowFlux;                        //Flux at the inflow boundary divided by 2PI

//Usefull position functions
extern double r_f(const Vector &x);                     //Function for r
extern double r_inv_f(const Vector &x);                 //Function for 1/r
extern void zero_f(const Vector &x, Vector &f);         //Function for 0 (vector)
extern void r_inv_hat_f(const Vector &x, Vector &f);    //Function for (1/r)*r^ (r^ unitary vector)
extern void rot_f(const Vector &x, DenseMatrix &f);     //Function for ( 0   1 )
                                                        //             (-1   0 )

//Physical properties (in T,S)
extern double FusionPoint(const double S);                              //Fusion temperature at a given salinity
extern double Phase(const double T, const double S);                    //Phase indicator (1 for liquid and 0 for solid)
extern double HeatInertia(const double T, const double S);              //Coefficient for the mass term in the temperature equation
extern double HeatDiffusivity(const double T, const double S);          //Coefficient for the diffusion term in the temperature equation
extern double SaltDiffusivity(const double T, const double S);          //Coefficient for the diffusion term in the salinity equation
extern double Impermeability(const double T, const double S);           //Inverse of the brinkman penalization permeability
extern double Density(const double T, const double S);                  //Relative density of the fluid
