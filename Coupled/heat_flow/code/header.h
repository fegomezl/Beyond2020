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
    double EpsilonEta;
    double c_l, c_s;
    double k_l, k_s;
    double L;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Vector &X);

        void SetParameters(const Vector &X, const Vector &rV);   //Update parameters from previous step

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);   //Sundials solver

        virtual ~Conduction_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> ess_tdof_list;

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator

        HypreParMatrix *M, *M_e, *M_0;
        HypreParMatrix *K_0;
        HypreParMatrix *T, *T_e;
        mutable HypreParVector Z;   

        //Solver objects
        HyprePCG M_solver;
        HyprePCG T_solver;
        HypreBoomerAMG M_prec;
        HypreBoomerAMG T_prec; 

        //Auxiliar grid functions
        ParGridFunction aux;
        ParGridFunction aux_C;
        ParGridFunction aux_K;
        ParGridFunction aux_L;
        ParGridFunction rv;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;

        ProductCoefficient coeff_rC;
        ProductCoefficient coeff_rK;

        VectorGridFunctionCoefficient coeff_rV;
        ScalarVectorProductCoefficient coeff_rCV;

        InnerProductCoefficient dHdT;
        InnerProductCoefficient dT_2;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const Vector &Theta);

        void SetParameters(const Vector &Theta);

        void Solve(Vector &W, Vector &Psi, Vector &V, Vector &rV);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> block_true_offsets;
        Array<int> ess_bdr_w, ess_bdr_psi;

        //System objects
        ParGridFunction psi;
        ParGridFunction w;
        ParGridFunction v;
        ParGridFunction rv;

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

        //Aditional variables
        ParGridFunction theta;
        ParGridFunction theta_eta;
        ParGridFunction psi_grad;
        ParGridFunction theta_dr;

        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient inv_R;
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

        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;
        HYPRE_Int size_v;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;

        HypreParVector *Theta;
        HypreParVector *W;
        HypreParVector *Psi;
        HypreParVector *rV;
        HypreParVector *V;

        //Operators
        Conduction_Operator *cond_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print objects
        ParaViewDataCollection *paraview_out;
};

extern double Rmin, Rmax, Zmin, Zmax;

extern double r_f(const Vector &x);
extern double inv_r(const Vector &x);
extern void zero_f(const Vector &x, Vector &f);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void r_inv_hat_f(const Vector &x, Vector &f);
