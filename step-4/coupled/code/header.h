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
    double epsilon_eta;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X);

        void SetParameters(const BlockVector &X);
        void UpdateVelocity(const ParGridFunction *v_flow);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

        virtual ~Conduction_Operator();
    protected:
        //Global parameters
        Config config;

        ParFiniteElementSpace &fespace;
        Array<int> block_true_offsets;

        mutable Array<int> newmann_bdr_theta, newmann_bdr_phi;
        mutable Array<int> robin_bdr_theta, robin_bdr_phi;

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

        ParGridFunction v;

        //Coefficients
        mutable FunctionCoefficient coeff_r;
        VectorFunctionCoefficient zero;
        MatrixFunctionCoefficient rot;

        VectorGridFunctionCoefficient coeff_rV;
        ScalarVectorProductCoefficient dt_coeff_rV;

        ProductCoefficient coeff_rCL;

        ProductCoefficient coeff_rK; 
        ProductCoefficient dt_coeff_rK;

        ProductCoefficient coeff_rD;
        ProductCoefficient dt_coeff_rD;

        ScalarVectorProductCoefficient coeff_rCLV;
        ScalarVectorProductCoefficient dt_coeff_rCLV;

        mutable FunctionCoefficient newmann_theta, newmann_phi;
        mutable ProductCoefficient r_newmann_theta, r_newmann_phi;
        mutable ProductCoefficient dt_r_newmann_theta, dt_r_newmann_phi;

        mutable FunctionCoefficient robin_h_theta, robin_h_phi;
        mutable FunctionCoefficient robin_ref_theta, robin_ref_phi;
        mutable ProductCoefficient neg_robin_h_theta, neg_robin_h_phi;
        mutable ProductCoefficient neg_robin_h_ref_theta, neg_robin_h_ref_phi;
        mutable ProductCoefficient neg_r_robin_h_theta, neg_r_robin_h_phi;
        mutable ProductCoefficient neg_r_robin_h_ref_theta, neg_r_robin_h_ref_phi;
        mutable ProductCoefficient neg_dt_r_robin_h_theta, neg_dt_r_robin_h_phi;
        mutable ProductCoefficient neg_dt_r_robin_h_ref_theta, neg_dt_r_robin_h_ref_phi;
};

class Flow_Operator{
  public:
    Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const BlockVector X);
    void Solve(const BlockVector X);
    void Update_T(const BlockVector X);
    ParGridFunction GetStream(){return *psi;}
    ParGridFunction GetVelocity(){return *v;}
    ParGridFunction GetVorticity(){return *w;}
    ~Flow_Operator();

  protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace;

       Array<int> block_true_offsets;

        //System objects
        ParLinearForm *f;
        ParLinearForm *g;
        ParBilinearForm *m;
        ParBilinearForm *d;
        ParMixedBilinearForm *c;
        Array<int> ess_bdr_psi;
        Array<int> ess_bdr_w;
        Array<int> ess_tdof_list_w;
        Array<int> ess_tdof_list_psi;

        //Solver objects
        ParGridFunction *psi;
        ParGridFunction *w;
        ParGridFunction *v;
        BlockVector Y;
        BlockVector B;
        HypreParMatrix *M;
        HypreParMatrix *D;
        HypreParMatrix *C;
        Array2D<HypreParMatrix*> hBlocks;
        HypreParMatrix *H;

       //Boundary conditions
       ParGridFunction *w_aux;
       ParGridFunction *psi_aux;
       ParGridFunction *v_aux;
       ParGridFunction *theta;
       ParGridFunction *phi;

       //Rotational coefficients
       FunctionCoefficient r;
       FunctionCoefficient r_invCoeff;
       VectorFunctionCoefficient r_hatCoeff;
       ScalarVectorProductCoefficient r_inv_hat;
       VectorFunctionCoefficient zero;
       //Properties coefficients
       GridFunctionCoefficient eta;
       ConstantCoefficient neg;
       //Dirichlet coefficients
       FunctionCoefficient w_coeff;
       VectorFunctionCoefficient w_grad;

       FunctionCoefficient psi_coeff;
       VectorFunctionCoefficient psi_grad;
       //Rotational coupled coefficients
       ScalarVectorProductCoefficient eta_r_inv_hat;

       ProductCoefficient neg_w;
       InnerProductCoefficient r_inv_hat_w_grad;
       ScalarVectorProductCoefficient neg_w_grad;

       InnerProductCoefficient r_inv_hat_psi_grad;
       InnerProductCoefficient eta_r_inv_hat_psi_grad;
       ScalarVectorProductCoefficient neg_psi_grad;
       ScalarVectorProductCoefficient eta_psi_grad;
       ScalarVectorProductCoefficient neg_eta_psi_grad;
       //Temperature coefficients
       DiscreteLinearOperator grad;
       MatrixFunctionCoefficient rot;
       GradientGridFunctionCoefficient gradpsi;
       MatrixVectorProductCoefficient rot_psi_grad;
       VectorGridFunctionCoefficient rV_aux;
       MatrixVectorProductCoefficient rV;

       //Solver objects
       SuperLUSolver *superlu;
       Operator *SLU_A;

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

        Array<int> block_true_offsets;

        //System objects
        ParGridFunction *theta;
        HypreParVector *Theta;
        ParGridFunction *phi;
        ParGridFunction *phase;

        BlockVector X;
        Conduction_Operator *oper_T;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        ParaViewDataCollection *paraview_out;

        //Flow_Operator objects
        Flow_Operator *flow_oper;
        ParGridFunction *x_psi;
        ParGridFunction *x_v;
        ParGridFunction *x_w;
};

extern double r_f(const Vector &x);
extern void rot_f(const Vector &x, DenseMatrix &f);
extern void zero_f(const Vector &x, Vector &f);

extern double T_fun(const double &salinity);

extern double delta_c_s_fun(const double &temperature, const double &salinity);

extern double Rmin, Rmax, Zmin, Zmax, height;
extern double epsilon_r;

extern double vel;
extern double entrance;

extern double theta_in, theta_out;
extern double phi_in, phi_out;
extern double rad;
