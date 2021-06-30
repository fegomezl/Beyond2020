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
    int order;
    int refinements;
    double dt_init;
    double t_final;
    int vis_steps_max;
    int ode_solver_type;
    double reltol;
    double abstol;
    double T_f;
    double c_s, c_l;
    double k_s, k_l;
    double L;
    double invDeltaT;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int attributes, Vector &X);

        void SetParameters(const Vector &X);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt); //Solver for implicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

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

        FunctionCoefficient r;

        GridFunctionCoefficient coeff_C;
        ProductCoefficient coeff_rC;

        GridFunctionCoefficient coeff_K;
        ProductCoefficient coeff_rK; ProductCoefficient dt_coeff_rK;

        GridFunctionCoefficient coeff_L;
        ProductCoefficient coeff_rL;
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
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;
};

extern double rf(const Vector &x);

extern double Rmin, Rmax, Zmin, Zmax;
