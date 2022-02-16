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

    double a_l, a_s;
    double sigma, kappa, eta;
};

class Transport_Operator : public TimeDependentOperator{        //Evolution of the conserved quantities (entalphy)
    public:
        Transport_Operator(Config config, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Vector &X);

        void SetParameters(const Vector &X, const Vector &r_Velocity);    //Update parameters from previous step

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);   //Sundials solver

        virtual ~Transport_Operator();
    protected:
        //Global parameters
        Config config;
 
        ParFiniteElementSpace &fespace_L2;
        ParFiniteElementSpace &fespace_ND;
        Array<int> ess_bdr;

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator
        ParLinearForm *b;    //RHS

        HypreParMatrix *M, *K ,*T;
        HypreParVector *B, *B_dt;
        mutable HypreParVector Z;   

        HypreBoomerAMG M_prec;
        HypreBoomerAMG T_prec;
        HyprePCG M_solver;
        GMRESSolver T_solver;

        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient coeff_zero;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1);

        void SetParameters(const Vector &X);

        void Solve(BlockVector &Y, Vector &Velocity, Vector &r_Velocity);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace_H1;
        ParFiniteElementSpace &fespace_L2;
        ParFiniteElementSpace &fespace_ND;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;

        //System objects
        ParGridFunction vorticity;
        ParGridFunction stream;

        //Solver objects
        HypreParVector Vorticity;
        HypreParVector Stream;

        HypreParMatrix *A00, *A00_e;
        HypreParMatrix *A11, *A11_e;
        HypreParMatrix *A01, *A01_e;
        HypreParMatrix *A10, *A10_e;

        HypreParVector B0;
        HypreParVector B1;
      
        //Rotational coefficients
        FunctionCoefficient coeff_r;
        FunctionCoefficient coeff_r_inv;
        VectorFunctionCoefficient coeff_r_inv_hat;
        MatrixFunctionCoefficient coeff_rot;

        //Construction of r_Velocity
        ParDiscreteLinearOperator gradient;
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

        FiniteElementCollection *fec_H1;
        FiniteElementCollection *fec_L2;
        FiniteElementCollection *fec_ND;

        ParFiniteElementSpace *fespace_H1;
        ParFiniteElementSpace *fespace_L2;
        ParFiniteElementSpace *fespace_ND;

        Array<int> block_offsets_H1;

        int dim;
        double h_min;
        int serial_refinements;

        HYPRE_Int size_H1;
        HYPRE_Int size_L2;
        HYPRE_Int size_ND;

        //System objects
        ParGridFunction *enthalphy;
        ParGridFunction *vorticity;
        ParGridFunction *stream;
        ParGridFunction *velocity;
        ParGridFunction *r_velocity;

        HypreParVector *X;          //Enthalphy
        BlockVector Y;              //Vorticity and Stream
        HypreParVector *Velocity;
        HypreParVector *r_Velocity;

        const IntegrationRule *irs[Geometry::NumGeom];

        //Operators
        Transport_Operator *transport_oper;
        Flow_Operator *flow_oper;

        //Solver objects
        ODESolver *ode_solver;
        ARKStepSolver *arkode;

        //Print objects
        ParaViewDataCollection *paraview_out;
};

//Simulation parameters
extern double R, Z;
extern double JumpScale;
extern double Epsilon;

//Rotational functions
extern double r_f(const Vector &x);
extern double r_inv_f(const Vector &x);
extern void zero_f(const Vector &x, Vector &f);
extern void r_inv_hat_f(const Vector &x, Vector &f);
extern void rot_f(const Vector &x, DenseMatrix &f);

//Relationship between variables
extern double TtoH(const double T);  //Map from temperature to enthaphy
extern double HtoT(const double H);  //Map from enthalphy to temperature
extern double Phase(const double H);    //Phase indicator (1 for liquid and 0 for solid)
extern double Density(const double H);  //Relative density of a corresponding state multiplied by g/v
