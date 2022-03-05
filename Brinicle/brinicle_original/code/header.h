#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

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

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);

        virtual ~Transport_Operator();
    protected:
        //All the variables with a 0 correspond to temperature
        //All the variables with a 1 correspond to salinity

        //Global parameters
        Config config;

        ParFiniteElementSpace &fespace_H1;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;
        Array<int> ess_bdr_0, ess_bdr_1;

        //Auxiliar grid functions
        ParGridFunction temperature, salinity, phase;
        ParGridFunction heat_inertia;
        ParGridFunction heat_diffusivity;
        ParGridFunction salt_diffusivity;
        ParGridFunction rvelocity;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient coeff_zero;

        //System objects
        HypreParMatrix *M0,   *M1;  
        HypreParMatrix *M0_e, *M1_e;
        HypreParMatrix *M0_o, *M1_o;
        HypreParMatrix *K0,   *K1;
        HypreParMatrix *T0,   *T1;
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
        //All the variables with a 0 correspond to vorticity
        //All the variables with a 1 correspond to stream

        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace_H1;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;
        Array<int> ess_bdr_0, ess_bdr_1;

        //System objects
        ParGridFunction velocity;
        ParGridFunction rvelocity;

        //Solver objects
        HypreParVector Vorticity, Stream;

        HypreParMatrix *A00,  *A00_e;
        HypreParMatrix *A01,  *A01_e;
        HypreParMatrix *A10,  *A10_e;
        HypreParMatrix *A11,  *A11_e;

        HypreParVector *B0, *B1;
        BlockVector B;

        //Additional variables
        ParGridFunction temperature;
        ParGridFunction temperature_dr;
        ParGridFunction salinity;
        ParGridFunction salinity_dr;
        ParGridFunction phase;
        ParGridFunction impermeability;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient coeff_r_inv;
        VectorFunctionCoefficient coeff_r_inv_hat;
        MatrixFunctionCoefficient coeff_rot;

        //Construction rV
        ParGridFunction stream;
        ParGridFunction stream_gradient;
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
        int vis_impressions;
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
        ParGridFunction *vorticity;
        ParGridFunction *stream;
        ParGridFunction *velocity;
        ParGridFunction *rvelocity;
        ParGridFunction *phase;

        BlockVector X;              //theta and phi
        BlockVector Y;              //w and psi
        HypreParVector *rVelocity;
        HypreParVector *Velocity;

        //Operators
        Transport_Operator *transport_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print parameters
        ParaViewDataCollection *paraview_out;
};

//Simulation parameters
extern double RMin, RMax, ZMin, ZMax;       //Size of the domain
extern double RIn, ZOut;                    //Size of the inflow and outflow
extern double Epsilon, EpsilonInv;          //Size of the indetermination window in heaviside functions 
                                            //(10̣^(-n) and 10^(n) respectively)

extern double InflowVelocity;                    //Velocity of the inflow
extern double InitialTemperature;                //Initial temperature of the domain
extern double InflowTemperature;                 //Temperature of the inflow
extern double InitialSalinity;                   //Initial salinity of the domain
extern double InflowSalinity;                    //Salinity of the inflow
extern double NucleationLength;                  //Lenght of the nucleation point
extern double NucleationHeight;                  //Height of the nucleation point
extern double NucleationTemperature;             //Temperature of the nucleation point
extern double NucleationSalinity;                //Salinity of the nucleation point
                                                 
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
