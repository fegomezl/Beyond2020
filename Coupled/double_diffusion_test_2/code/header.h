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

    double c_l, c_s;
    double k_l, k_s;
    double d_l, d_s;
    double L_l, L_s;
    double invDeltaT;
    double EpsilonT;
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
        Array<int> ess_tdof_list;

        //System objects
        SuperLURowLocMatrix *M, *T;
        HypreParMatrix *M_0, *K_0, *T_e;

        mutable BlockVector Z; 
        BlockVector dZ, T_d;

        //Solver objects
        SuperLUSolver *M_solver, *T_solver;

        //Auxiliar grid functions
        ParGridFunction phi, theta, phase;
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

//Simulation parameters
extern double Rmin, Rmax, Zmin, Zmax;

//Rotational functions
extern double r_f(const Vector &x);
extern void zero_f(const Vector &x, Vector &f);

//Fusion temperature dependent of salinity
extern double T_fun(const double &salinity);
