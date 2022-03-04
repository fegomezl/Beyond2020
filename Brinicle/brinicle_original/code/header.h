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
        ParGridFunction aux_C;
        ParGridFunction aux_K;
        ParGridFunction aux_D;
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

        HypreParVector B0, B1;
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
        Array<int> ess_tdof_w, ess_tdof_psi;
        Array<int> ess_bdr_w, ess_bdr_psi;
        Array<int> bdr_psi_in, bdr_psi_out;
        Array<int> bdr_psi_closed_down, bdr_psi_closed_up;

        //System objects
        ParGridFunction v;
        ParGridFunction rv;

        //Solver objects
        HypreParVector W, Psi;

        HypreParMatrix *M,  *M_e;
        HypreParMatrix *D,  *D_e;
        HypreParMatrix *C,  *C_e;
        HypreParMatrix *Ct, *Ct_e;

        HypreParVector F, G;
        BlockVector B;

        //Additional variables
        ParGridFunction theta;
        ParGridFunction theta_dr;
        ParGridFunction phi;
        ParGridFunction phi_dr;
        ParGridFunction phase;
        ParGridFunction eta;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient inv_R;
        VectorFunctionCoefficient r_inv_hat;
        MatrixFunctionCoefficient rot;

        //Construction rV
        ParGridFunction psi;
        ParGridFunction psi_grad;
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
        int vis_impressions;
        double total_time;

        //Mesh objects
        ParMesh *pmesh;

        FiniteElementCollection *fec;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace;
        ParFiniteElementSpace *fespace_v;

        Array<int> block_true_offsets;
        
        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;
        HYPRE_Int size_v;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *phi;
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;
        ParGridFunction *rv;
        ParGridFunction *phase;

        BlockVector X;              //theta and phi
        BlockVector Z;              //w and psi
        HypreParVector *rV;
        HypreParVector *V;

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

//Usefull position functions
extern double r_f(const Vector &x);                     //Function for r
extern double r_inv_f(const Vector &x);                 //Function for 1/r
extern void zero_f(const Vector &x, Vector &f);         //Function for 0 (vector)
extern void r_inv_hat_f(const Vector &x, Vector &f);    //Function for (1/r)*r^ (r^ unitary vector)
extern void rot_f(const Vector &x, DenseMatrix &f);     //Function for ( 0   1 )
                                                        //             (-1   0 )
//Physical Properties
extern double T_fun(const double &salinity);
extern double delta_c_s_fun(const double &temperature, const double &salinity);
extern double delta_k_s_fun(const double &temperature, const double &salinity);
extern double delta_l_s_fun(const double &temperature, const double &salinity);
extern double delta_rho_t_fun(const double &temperature, const double &salinity);
extern double delta_rho_p_fun(const double &temperature, const double &salinity);

//Brinicle conditions
extern double InflowVelocity;           //Velocity of the inflow
extern double InflowFlux;               //Flux of the given inflow velocity field
extern double InitialTemperature;       //Initial temperature of the domain
extern double InflowTemperature;        //Temperature of the inflow
extern double InitialSalinity;          //Initial salinity of the domain
extern double InflowSalinity;           //Salinity of the inflow
extern double NucleationLength;         //Lenght of the nucleation point
extern double NucleationHeight;         //Height of the nucleation point
extern double NucleationTemperature;    //Temperature of the nucleation point
extern double NucleationSalinity;       //Salinity of the nucleation point
