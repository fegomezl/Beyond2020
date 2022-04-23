#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Constants{
    //Properties for the system 

    //Fusion point coeficients (T_f = -(aS + BS^3))
    double FusionPoint_a = 6.04E-1; 
    double FusionPoint_b = 5.81E-4;

    //Heat inertia (Specific heat capacity over latent heat c/L) (^C)
    double m_l = 1.17E-2;    //Liquid
    double m_s = 6.57E-3;    //Solid

    //Thermal diffusivity (Thermal conductivity over volumetric latent heat k/(rho*L)) (cm^2/s)
    double d_temperature_l = 1.64E-5;   //Liquid
    double d_temperature_s = 7.45E-5;   //Solid

    //Salt diffusivity (cm^2/s)
    double d_salinity_l = 1.67E-5;      //Liquid
    double d_salinity_s = 0;            //Solid

    //Density coeficients (\Delta\rho)
    double Buoyancy_k = 1.44E4;          //g/nu (1/(cm*s)) 
    double Buoyancy_a0 = 8.25E-3;        //T^(0)*S^(1)
    double Buoyancy_a1 = -4.11E-5;       //T^(1)*S^(1)
    double Buoyancy_a2 = 7.73E-7;        //T^(2)*S^(1)
    double Buoyancy_a3 = -8.32E-9;       //T^(3)*S^(1)
    double Buoyancy_a4 = 5.52E-11;       //T^(4)*S^(1)
    double Buoyancy_b0 = -1.83E-4;       //T^(0)*S^(1.5)
    double Buoyancy_b1 = 3.33E-6;        //T^(1)*S^(1.5)
    double Buoyancy_b2 = -5.60E-8;       //T^(2)*S^(1.5)
    double Buoyancy_c0 = 4.89E-5;        //T^(0)*S^(2)
};

struct Config{
    Config(int pid, int nproc);

    bool master;
    int nproc;
    int pid;

    double dt_init;
    double t_final;
    int vis_steps_max;
    bool rescale;

    int refinements;
    int order;
    double reltol_conduction;
    double abstol_conduction;
    int iter_conduction;
    double reltol_sundials;
    double abstol_sundials;

    bool restart;
    double t_init;
};

class Transport_Operator : public TimeDependentOperator{
    public:
        Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X);

        void SetParameters(const BlockVector &X, const Vector &rVelocity);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);

        virtual ~Transport_Operator();
    protected:
        //All 0-variables are related to temperature
        //All 1-variables are related to salinity

        //Global parameters
        Config config;

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

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X);

        void SetParameters(const BlockVector &X);

        void Solve(BlockVector &Y, Vector &Velocity, Vector &rVelocity);

        ~Flow_Operator();
    protected:
        //All 0-variables are related to vorticity
        //All 1-variables are related to stream

        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace_H1;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;
        Array<int> ess_bdr_0, ess_bdr_1;
        Array<int> ess_bdr_in, ess_bdr_out;
        Array<int> ess_bdr_closed_down, ess_bdr_closed_up;

        //System objects
        ParGridFunction vorticity_boundary;
        ParGridFunction stream_boundary;
        ParGridFunction velocity;
        ParGridFunction rvelocity;

        //Solver objects
        BlockVector B;
        HypreParVector *B0, *B1;

        HypreParMatrix *A00;
        HypreParMatrix *A01;
        HypreParMatrix *A10;
        HypreParMatrix *A11;

        //Additional variables
        ParGridFunction temperature;
        ParGridFunction temperature_dr;
        ParGridFunction salinity;
        ParGridFunction salinity_dr;
        ParGridFunction phase;
        ParGridFunction impermeability;
        ParGridFunction stream;
        ParGridFunction stream_gradient;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient coeff_r_inv;
        VectorFunctionCoefficient coeff_r_inv_hat;
        MatrixFunctionCoefficient coeff_rot;

        //Boundary coefficients
        FunctionCoefficient coeff_vorticity;
        FunctionCoefficient coeff_stream;
        FunctionCoefficient coeff_stream_in;
        FunctionCoefficient coeff_stream_out;
        ConstantCoefficient coeff_stream_closed_down;
        ConstantCoefficient coeff_stream_closed_up;

        //Construction rV
        ParDiscreteLinearOperator gradient;
};

class Artic_sea{
    public:
        Artic_sea(Config config);
        void run(const char *mesh_file);
        ~Artic_sea();
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void time_step();
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

        //Mesh objects
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

        //Operators
        Transport_Operator *transport_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print parameters
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
extern double TemperatureMax, TemperatureMin;    //Limits of the temperature scale
extern double SalinityMax, SalinityMin;          //Limits of the salinity scale

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
extern double HeatInertia(const double T, const double S);              //Heat capacity over latent heat
extern double HeatDiffusivity(const double T, const double S);          //Heat conduction over latent heat
extern double SaltDiffusivity(const double T, const double S);          //Diffusion coefficient for the mass equation
extern double Impermeability(const double T, const double S);           //Inverse of the brinkman penalization permeability
extern double ExpansivityTemperature(const double T, const double S);   //Expansivity coefficient for the temperature gradient
extern double ExpansivitySalinity(const double T, const double S);      //Expansivity coefficient for the salinity gradient

//Bounding variables
extern double T_bounded(const double T);        //Bounding of temperature according to expectations
extern double S_bounded(const double S);        //Bounding of salinity according to expectations
