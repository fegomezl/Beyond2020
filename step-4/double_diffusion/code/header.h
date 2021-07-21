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
    double D_l, D_s;
    double L;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X);

        void SetParameters(const BlockVector &X);
        void UpdateVelocity(const Vector &psi);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

        virtual ~Conduction_Operator();
    protected:
        //Global parameters
        Config config;
        double Scaled_dt;

        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> block_true_offsets;

        mutable Array<int> newmann_bdr_theta, newmann_bdr_phi;

        Array<int> ess_tdof_list_theta, ess_tdof_list_phi;

        //System objects
        ParBilinearForm *m_theta, *m_phi;        //Mass operators
        ParBilinearForm *k_theta, *k_phi;        //Difussion operators
        ParBilinearForm *t_theta, *t_phi;        //m + dt*k

        HypreParMatrix M_theta, M_phi;
        HypreParMatrix T_theta, T_phi;   

        //Solver objects
        CGSolver M_theta_solver, M_phi_solver;
        CGSolver T_theta_solver, T_phi_solver;
        HypreSmoother M_theta_prec, M_phi_prec;
        HypreSmoother T_theta_prec, T_phi_prec;

        //Auxiliar grid functions
        ParGridFunction aux_phi, aux_theta;
        ParGridFunction aux_C;
        ParGridFunction aux_K;
        ParGridFunction aux_D;

        ParGridFunction psi;

        //Coefficients
        mutable FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;
        MatrixFunctionCoefficient rot;

        GradientGridFunctionCoefficient gradpsi;
        MatrixVectorProductCoefficient coeff_rV;
        ScalarVectorProductCoefficient dt_coeff_rV;

        ProductCoefficient coeff_rCL;

        ProductCoefficient coeff_rK; 
        ProductCoefficient dt_coeff_rK;

        ProductCoefficient coeff_rD; 
        ProductCoefficient dt_coeff_rD;

        ScalarVectorProductCoefficient coeff_rCLV;
        ScalarVectorProductCoefficient dt_coeff_rCLV;

        mutable FunctionCoefficient newmann_theta, newmann_phi;
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

        Array<int> block_true_offsets;

        //System objects
        ParGridFunction *theta;
        ParGridFunction *phi;
        ParGridFunction *phase;

        BlockVector X;
        Conduction_Operator *oper_T;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;
};

extern double r_f(const Vector &x);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void zero_f(const Vector &x, Vector &f);

extern double T_fun(const double &salinity);

extern double delta_c_s_fun(const double &temperature, const double &salinity);

extern double Rmin, Rmax, Zmin, Zmax;

extern double mid;
