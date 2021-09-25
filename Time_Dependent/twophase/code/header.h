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
    double L;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Vector &X);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);   //Sundials solver

        void SetParameters(const Vector &X);    //Update parameters from previous step

        virtual ~Conduction_Operator();

 protected:
        //Global parameters
        Config config;

        ParFiniteElementSpace &fespace;
        Array<int> ess_tdof_list;

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator
        ParBilinearForm *t;  //m + dt*k

        HypreParMatrix M;
        HypreParMatrix K;
        HypreParMatrix T;   

        CGSolver M_solver;
        CGSolver T_solver;
        HypreSmoother M_prec;
        HypreSmoother T_prec; 

        ParGridFunction aux;
        ParGridFunction aux_C;
        ParGridFunction aux_K;
        ParGridFunction aux_L;

        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;

        ProductCoefficient coeff_rC;
        ProductCoefficient coeff_rK; ProductCoefficient dt_coeff_rK;

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

        int iteration;
        double t;
        double dt;
        bool last;
        int vis_iteration;
        int vis_steps;
        int vis_impressions;

        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;

        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParGridFunction *x;
        HypreParVector *X;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;
};

extern double initial(const Vector &x);

extern double Rmin, Rmax, Zmin, Zmax;
