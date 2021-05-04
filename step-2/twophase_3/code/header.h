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

};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr, Array<int> nbc_marker);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt,
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &b,
                                     int j_update, int *j_status, double scaled_dt);

	      virtual int SUNImplicitSolve(const Vector &b, Vector &X,
                                     double tol);

        void SetParameters(const Vector &X);

        virtual ~Conduction_Operator();

 protected:
        ParFiniteElementSpace &fespace;
        Array<int> ess_tdof_list;

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator
        ParLinearForm *f;    //RHS

        HypreParMatrix M;
        HypreParMatrix K;
        HypreParMatrix *T;    //T = M + dt K
        HypreParVector F;

        CGSolver M_solver;
        CGSolver T_solver;
        HypreSmoother M_prec;
        HypreSmoother T_prec;

        FunctionCoefficient r;
        VectorFunctionCoefficient r_hat;

        mutable HypreParVector z;
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
        Array<int> ess_bdr;
        Array<int> nbc_marker;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;

};

extern double T_f;
extern double Rmin, Rmax, Zmin, Zmax;

extern double alpha_l; //Liquid thermal conduction
extern double alpha_s; //Solid thermal conduction
