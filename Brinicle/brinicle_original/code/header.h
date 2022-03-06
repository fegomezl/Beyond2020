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

    double T_f;
    double c_l, c_s;
    double k_l, k_s;
    double d_l, d_s;
    double L_l, L_s;

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
        Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X);

        void SetParameters(const BlockVector &X);

        void Solve(BlockVector &Z, Vector &V, Vector &rV);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> block_true_offsets;
        Array<int> ess_bdr_w, ess_bdr_psi;
        Array<int> bdr_psi_in, bdr_psi_out;
        Array<int> bdr_psi_closed_down, bdr_psi_closed_up;

        //System objects
        ParGridFunction psi;
        ParGridFunction w;
        ParGridFunction v;
        ParGridFunction rv;

        ParLinearForm *f;
        ParLinearForm *g;
        ParBilinearForm *m;
        ParBilinearForm *d;
        ParMixedBilinearForm *c;
        ParMixedBilinearForm *ct;

        //Solver objects
        BlockVector Y;
        BlockVector B;

        HypreParMatrix *M;
        HypreParMatrix *D;
        HypreParMatrix *C;
        HypreParMatrix *Ct;

        //Additional variables
        ParGridFunction theta;
        ParGridFunction theta_dr;
        ParGridFunction phi;
        ParGridFunction phi_dr;
        ParGridFunction phase;
        ParGridFunction eta;
        ParGridFunction psi_grad;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient inv_R;
        VectorFunctionCoefficient r_inv_hat;
        MatrixFunctionCoefficient rot;

        //Boundary coefficients
        FunctionCoefficient w_coeff;
        FunctionCoefficient psi_coeff;
        FunctionCoefficient psi_in;
        FunctionCoefficient psi_out;
        ConstantCoefficient closed_down;
        ConstantCoefficient closed_up;

        //Construction rV
        ParDiscreteLinearOperator grad;

        VectorGridFunctionCoefficient Psi_grad;
        MatrixVectorProductCoefficient rot_Psi_grad;
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
