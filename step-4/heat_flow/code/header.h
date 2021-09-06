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
    double epsilon_eta;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Vector &X);

        void SetParameters(const Vector &X, const Vector &V);

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
        ParGridFunction v;

        FunctionCoefficient r;
        VectorFunctionCoefficient zero;

        ProductCoefficient coeff_rCL;

        ProductCoefficient coeff_rK; 
        ProductCoefficient dt_coeff_rK;

        //GradientGridFunctionCoefficient rV;
        VectorGridFunctionCoefficient rV;
        ScalarVectorProductCoefficient coeff_rCLV;
        ScalarVectorProductCoefficient dt_coeff_rCLV;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const Vector &Theta);

        void SetParameters(const Vector &Theta);

        void Solve(Vector &W, Vector &Psi, Vector &V);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        ParGridFunction psi;
        ParGridFunction w;
        ParGridFunction v;

        //Mesh objects
        ParFiniteElementSpace &fespace;

        Array<int> block_true_offsets;

        Array<int> ess_bdr_w;
        Array<int> ess_bdr_psi;

        //System objects
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
        HypreParMatrix *Ct;

        Array2D<HypreParMatrix*> hBlocks;
        Array2D<double> blockCoeff; 
        HypreParMatrix *H;
       
        SuperLUSolver superlu;
        SuperLURowLocMatrix *SLU_A;

        //Aditional variables
        ParGridFunction theta;
        ParGridFunction theta_eta;
        ParGridFunction psi_grad;
        ParGridFunction theta_dr;

        //Rotational coefficients
        FunctionCoefficient r;
        VectorFunctionCoefficient r_inv_hat;
        MatrixFunctionCoefficient rot;

        //Boundary Coefficients
        FunctionCoefficient w_coeff; 
        FunctionCoefficient psi_coeff; 

        //Construction rV
        DiscreteLinearOperator grad;

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
        HYPRE_Int size_v;

        ParMesh *pmesh;

        FiniteElementCollection *fec;
        FiniteElementCollection *fec_v;

        ParFiniteElementSpace *fespace;
        ParFiniteElementSpace *fespace_v;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;

        HypreParVector *Theta;
        HypreParVector *W;
        HypreParVector *Psi;
        HypreParVector *V;

        //Operators
        Conduction_Operator *cond_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;
};

extern void zero_f(const Vector &x, Vector &f);
extern double r_f(const Vector &x);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void r_inv_hat_f(const Vector &x, Vector &f);

extern double Rmin, Rmax, Zmin, Zmax;
extern double epsilon_r;
