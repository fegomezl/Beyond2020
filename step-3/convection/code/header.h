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
        Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt,
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
	      virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

        void SetParameters(const Vector &X);

        virtual ~Conduction_Operator();

 protected:
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
        double h_min;

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
        FunctionCoefficient initial_f;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;

        //Flow_Operator objects
        FiniteElementCollection *fec_psi;
        FiniteElementCollection *fec_w;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace_psi;
        ParFiniteElementSpace *fespace_w;
        ParFiniteElementSpace *fespace_v;


};

class Flow_Operator{
  public:
        Flow_Operator(ParFiniteElementSpace &fespace_psi, ParFiniteElementSpace fespace_w,ParFiniteElementSpace fespace_v);
        ~Flow_Operator();

  protected:

        //Mesh objects

        ParFiniteElementSpace &fespace_psi;
        ParFiniteElementSpace &fespace_w;
        ParFiniteElementSpace &fespace_v;

        Array<int> block_offsets;
        Array<int> block_true_offsets;

        //System objects
        BlockVector x;
        BlockVector b;
        ParLinearForm *f;
        ParLinearForm *g;
        ParBilinearForm *m;
        ParBilinearForm *d;
        ParMixedBilinearForm *c;
        ParMixedBilinearForm *ct;

        //Solver objects
        BlockVector X;
        BlockVector B;
        HypreParMatrix *M;
        HypreParMatrix *D;
        HypreParMatrix *C;
        //TransposeOperator *Ct;
        HypreParMatrix *Ct;
        BlockOperator *A;
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;


};

extern double initial(const Vector &x);

extern double T_f;
extern double Rmin, Rmax, Zmin, Zmax;

extern double c_s, c_l;
extern double k_s, k_l;
extern double L;

extern double DeltaT;
