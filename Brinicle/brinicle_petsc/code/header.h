#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Config{
    Config(int pid, int nproc);

    bool master;
    int nproc;
    int pid;

    double dt_init;
    double t_final;
    int vis_steps_max;
    bool rescale;

    int refinements;
    int order;
    double reltol_conduction;
    double abstol_conduction;
    int iter_conduction;
    double reltol_sundials;
    double abstol_sundials;

    double invDeltaT;
    double EpsilonT;
    double EpsilonEta;
    double c_l, c_s;
    double k_l, k_s;
    double d_l, d_s;
    double L_l, L_s;

    bool restart;
    double t_init;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X);

        void SetParameters(const BlockVector &X, const Vector &rV);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

        virtual ~Conduction_Operator();
    protected:
        //Global parameters
        Config config;

        ParFiniteElementSpace &fespace, &fespace_v;
        Array<int> block_true_offsets;
        Array<int> ess_tdof_theta, ess_tdof_phi;

        //System objects
        HypreParMatrix *M_theta, *M_e_theta, *M_0_theta,    *M_phi, *M_e_phi, *M_0_phi;
        HypreParMatrix *K_0_theta,                          *K_0_phi;
        HypreParMatrix *T_theta, *T_e_theta,                *T_phi, *T_e_phi;

        mutable HypreParVector Z_theta, Z_phi;

        //Solver objects
        HyprePCG M_theta_solver, M_phi_solver;
        HyprePCG T_theta_solver, T_phi_solver;
        HypreBoomerAMG M_theta_prec, M_phi_prec;
        HypreBoomerAMG T_theta_prec, T_phi_prec;

        //Coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets);

        void SetParameters(const BlockVector &X);

        void Solve(BlockVector &Z, Vector &V, Vector &rV);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace, &fespace_v;
        Array<int> block_true_offsets;
        Array<int> ess_tdof_w, ess_tdof_psi;

        //System objects
        ParGridFunction psi;
        ParGridFunction w;

        //Solver objects
        HypreParVector W;
        HypreParVector Psi;
        HypreParVector B_w;
        HypreParVector B_psi;

        HypreParMatrix *M,  *M_e;
        HypreParMatrix *D,  *D_e;
        HypreParMatrix *C,  *C_e;
        HypreParMatrix *Ct, *Ct_e;

        PetscLinearSolver solver;
        PetscFieldSplitSolver prec;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient r_inv_hat;
        MatrixFunctionCoefficient rot;
        FunctionCoefficient inv_R;

        //Construction rV
        ParDiscreteLinearOperator grad;
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

        Array<int> block_true_offsets;
        
        int dim;
        double h_min;
        int serial_refinements;
        HYPRE_Int size;
        HYPRE_Int size_v;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *phi;
        ParGridFunction *w;
        ParGridFunction *psi;
        ParGridFunction *v;
        ParGridFunction *rv;
        ParGridFunction *phase;

        BlockVector X;              //theta and phi
        BlockVector Z;              //w and psi
        HypreParVector *rV;
        HypreParVector *V;

        const IntegrationRule *irs[Geometry::NumGeom];

        //Operators
        Conduction_Operator *cond_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print parameters
        ParaViewDataCollection *paraview_out;
};

//Simulation parameters
extern double Rmin, Rmax, Zmin, Zmax;
extern double R_in, Z_out;

//Rotational functions
extern double r_f(const Vector &x);
extern double inv_r(const Vector &x);
extern void zero_f(const Vector &x, Vector &f);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void r_inv_hat_f(const Vector &x, Vector &f);

//Fusion temperature dependent of salinity
extern double T_fun(const double &salinity);

//Parameters of buoyancy
extern double delta_rho_t_fun(const double &temperature, const double &salinity);
extern double delta_rho_p_fun(const double &temperature, const double &salinity);

//Brinicle conditions
extern double Vel, Q;
extern double theta_in, theta_out;
extern double phi_in, phi_out;
extern double n_l, n_h;
extern double theta_n, phi_n;
extern double c_l;