#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Config{
    Config(bool master, int nproc);

    bool master;
    int nproc;

    double dt_init;
    double t_final;
    int vis_steps_max;

    int refinements;
    int order;
    double reltol_conduction;
    double abstol_conduction;
    int iter_conduction;
    double reltol_sundials;
    double abstol_sundials;

    double T_f;
    double invDeltaT;
    double EpsilonT;
    double c_l, c_s;
    double k_l, k_s;
    double D_l, D_s;
    double L;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X);

        void SetParameters(const BlockVector &X);       //Update parameters from previous step

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);

        virtual ~Conduction_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> block_true_offsets;
        Array<int> ess_tdof_theta, ess_tdof_phi;

        Array<int> robin_bdr_theta, robin_bdr_phi;

        //System objects
        ParBilinearForm *m_theta, *m_phi;        //Mass operators
        ParBilinearForm *k_theta, *k_phi;        //Difussion operators

        HypreParMatrix *M_theta, *M_e_theta, *M_0_theta,    *M_phi, *M_e_phi, *M_0_phi;
        HypreParMatrix *K_0_theta,                          *K_0_phi;
        HypreParMatrix *T_theta, *T_e_theta,                *T_phi, *T_e_phi;

        HypreParVector F_theta, F_phi;
        HypreParVector dt_F_theta, dt_F_phi;
        mutable HypreParVector Z_theta, Z_phi;

        //Solver objects
        HyprePCG M_theta_solver, M_phi_solver;
        HyprePCG T_theta_solver, T_phi_solver;
        HypreBoomerAMG M_theta_prec, M_phi_prec;
        HypreBoomerAMG T_theta_prec, T_phi_prec;

        //Auxiliar grid functions
        ParGridFunction aux_phi, aux_theta;
        ParGridFunction aux_C;
        ParGridFunction aux_K;
        ParGridFunction aux_D;
        ParGridFunction aux_L;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;

        ProductCoefficient coeff_rC;
        ProductCoefficient coeff_rK; 
        ProductCoefficient coeff_rD; 
        
        InnerProductCoefficient dHdT;
        InnerProductCoefficient dT_2;

        FunctionCoefficient robin_h_theta, robin_h_phi;
        ProductCoefficient r_robin_h_theta, r_robin_h_phi;
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
        ParFiniteElementSpace *fespace;
        Array<int> block_true_offsets;

        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *phi;
        ParGridFunction *phase;

        BlockVector X;
        Conduction_Operator *cond_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print objects
        ParaViewDataCollection *paraview_out;
};

extern double r_f(const Vector &x);
extern void zero_f(const Vector &x, Vector &f);

extern double T_fun(const double &salinity);

extern double delta_c_s_fun(const double &temperature, const double &salinity);

extern double Rmin, Rmax, Zmin, Zmax;
