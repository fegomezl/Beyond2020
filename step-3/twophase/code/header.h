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
    int ode_solver_type;
    double reltol_conduction;
    double abstol_conduction;
    int iter_conduction;
    double reltol_sundials;
    double abstol_sundials;

    double T_f;
    double invDeltaT;
    double c_l, c_s;
    double k_l, k_s;
    double L;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Vector &X);

        void SetParameters(const Vector &X);
        void UpdateVelocity(const HypreParVector &psi);

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
        HypreParMatrix T;   

        CGSolver M_solver;
        CGSolver T_solver;
        HypreSmoother M_prec;
        HypreSmoother T_prec;

        ParGridFunction aux;
        ParGridFunction aux_C;
        ParGridFunction aux_K;

        ParGridFunction psi;

        FunctionCoefficient r;
        VectorFunctionCoefficient zero;

        ProductCoefficient coeff_rCL;

        ProductCoefficient coeff_rK; 
        ProductCoefficient dt_coeff_rK;

        //GradientGridFunctionCoefficient rV;
        MatrixFunctionCoefficient rot;
        GradientGridFunctionCoefficient gradpsi;
        MatrixVectorProductCoefficient rV;
        ScalarVectorProductCoefficient coeff_rCLV;
        ScalarVectorProductCoefficient dt_coeff_rCLV;
};

class Flow_Operator{
  public:
    Flow_Operator(Config config, ParFiniteElementSpace &fespace, int attributes, const ParGridFunction *x_T);
    void Solve(Config config, HypreParVector *X_Psi, ParGridFunction *x_psi, const ParGridFunction *x_T);
    void Update_T(Config config, const ParGridFunction *x_T);
    ParGridFunction *psi;
    ~Flow_Operator();

  protected:

        //Mesh objects
        ParFiniteElementSpace &fespace;

        Array<int> block_offsets;
        Array<int> block_true_offsets;

        //System objects
        BlockVector y;
        BlockVector b;
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

        //TransposeOperator *Ct;
        HypreParMatrix *Ct;
        BlockOperator *A;

        ParGridFunction *w;

        ConstantCoefficient bg;

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
        ParGridFunction *x_T;
        HypreParVector *X_T;
        Conduction_Operator *oper_T;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;

        //Flow_Operator objects
        Flow_Operator *flow_oper;
        ParGridFunction *x_psi;
        HypreParVector *X_Psi;
};

extern double r_f(const Vector &x);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void zero_f(const Vector &x, Vector &f);

extern double Rmin, Rmax, Zmin, Zmax;
