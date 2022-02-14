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

    double a_l, a_s;
    double sigma, kappa, eta;
};

class Transport_Operator : public TimeDependentOperator{        //Evolution of the conserved quantities (entalphy)
    public:
        Transport_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Vector &X);

        void SetParameters(const Vector &X);    //Update parameters from previous step

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);   //Sundials solver

        virtual ~Transport_Operator();
    protected:
        //Global parameters
        Config config;
 
        ParFiniteElementSpace &fespace;
        Array<int> ess_bdr, ess_bdr_l, ess_bdr_s;

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator
        ParLinearForm *b;    //RHS

        HypreParMatrix *M, *K ,*T;
        HypreParVector *B, *B_dt;
        mutable HypreParVector Z;   

        HypreBoomerAMG M_prec;
        HypreBoomerAMG T_prec;
        HyprePCG M_solver;
        GMRESSolver T_solver;

        ParGridFunction aux;

        VectorFunctionCoefficient zero;
        FunctionCoefficient coeff_r;
        ProductCoefficient coeff_rA, coeff_rA_l, coeff_rA_s;
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

        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;

        //System objects
        ParGridFunction *x;
        HypreParVector *X;
        Transport_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print objects
        ParaViewDataCollection *paraview_out;
};

extern double R, Z;
extern double JumpScale;

extern double H_l, H_s;

extern double H(const double T);  //Map from temperature to enthaphy
extern double T(const double H);  //Map from enthalphy to temperature
