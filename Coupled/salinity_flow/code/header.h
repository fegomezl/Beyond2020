#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Constants{
    //Explicacion 

    //Forces properties
    double v = 408;             //Kinematic viscosity \nu (mm^2/min)
    double g = 35280000;        //Gravity (mm/min^2)

    //Densities (kg/m3)
    double rho_0 = 999.8;       //Density of pure water 
    double rho_l = 1028.3;      //Density of sea water
    double rho_s = 934.4;       //Density of sea ice 

    //Energy properties
    double c_l = 3.77;              //Specific heat capacity of sea water (kJ/(kg*°C))
    double c_s = 2.12;              //Specific heat capacity of sea ice (kJ/(kg*°C))
    double k_l = 0.5433;            //Thermal diffusivity of sea water (W/(m*°C))
    double k_s = 2.2467;            //Thermal diffusivity of sea ice (W/(m*°C))
    double L = 322.7;               //Latent Heat (kJ/kg)

    //Fusion point coeficients (°C)
    double FusionPoint_a = 0.6037;       
    double FusionPoint_b = 0.00058123;   

    //Heat inertia (Specific heat capacity over latent heat c/L) (Dimentionless)
    double m_l = c_l/L;    //Liquid
    double m_s = c_s/L;    //Solid

    //Thermal diffusivity (Thermal conductivity over volumetric latent heat k/(rho*L)) (mm^2/min)
    double d_temperature_l = k_l/(rho_l*L);      //Liquid
    double d_temperature_s = k_s/(rho_s/L);    //Solid

    //Salt diffusivity (mm^2/min)
    double d_salinity_l = 0.1;                 //Liquid
    double d_salinity_s = 0.000000001;         //Solid

    //Density coeficients (\Delta\rho)
    double Buoyancy_k = g/(v*rho_0); 
    double Buoyancy_a0 = 8.246121;              //T^(0)*S^(1)
    double Buoyancy_a1 = -0.04109344;           //T^(1)*S^(1)
    double Buoyancy_a2 = 0.0007731839;          //T^(2)*S^(1)
    double Buoyancy_a3 = -0.000008323176;       //T^(3)*S^(1)
    double Buoyancy_a4 = 0.00000005516023;      //T^(4)*S^(1)
    double Buoyancy_b0 = -0.182739;             //T^(0)*S^(1.5)
    double Buoyancy_b1 = 0.003327135;           //T^(1)*S^(1.5)
    double Buoyancy_b2 = -0.00005600674;        //T^(2)*S^(1.5)
    double Buoyancy_c0 = 0.04886169;            //T^(0)*S^(2)
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
        Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, Vector &Temperature, Vector &Salinity);

        void SetParameters(const Vector &Temperature, const Vector &Salinity, const Vector &rVelocity);

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
        Array<int> ess_tdof;
        Array<int> ess_bdr;

        //Auxiliar grid functions
        ParGridFunction temperature, salinity;
        ParGridFunction rvelocity;
        ParGridFunction salt_diffusivity;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient coeff_zero;

        ProductCoefficient coeff_rD; 

        VectorGridFunctionCoefficient coeff_rV;

        //System objects
        HypreParMatrix *M;
        HypreParMatrix *M_e; 
        HypreParMatrix *M_o;
        HypreParMatrix *K;
        HypreParMatrix *T;
        HypreParMatrix *T_e;

        HypreParVector *B;
        HypreParVector B_dt;
        mutable HypreParVector Z;

        //Solver objects
        HyprePCG M_solver;
        HyprePCG T_solver;
        HypreBoomerAMG M_prec;
        HypreBoomerAMG T_prec;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1);

        void SetParameters(const Vector &Temperature, const Vector &Salinity);

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

        HypreParVector *Temperature;
        HypreParVector *Salinity;
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
