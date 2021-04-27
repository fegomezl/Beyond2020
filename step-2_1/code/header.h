#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Config{
    //Constructor
    Config(bool master, int nproc);
    //Passing parameters
    bool master;
    int nproc;
    int order;
    int refinements;
    double dt_init;
    double t_init;
    double t_final;
    int vis_steps;
    int ode_solver_type;
    double reltol;
    double abstol;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(ParFiniteElementSpace &fespace, const HypreParVector &X, Array<int> ess_bdr, double t_init);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt, 
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &b, 
                                     int j_update, int *j_status, double scaled_dt);
	      virtual int SUNImplicitSolve(const Vector &b, Vector &X, 
                                     double tol);

        void SetParameters(const HypreParVector &X, Array<int> ess_bdr);                        //Update the bilinear forms

        virtual ~Conduction_Operator();
    protected:
        //Mesh objects
        ParFiniteElementSpace &fespace;
        Array<int> ess_tdof_list; 

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator
        ParLinearForm *f;    //Bounary term

        //Solver objects
        HypreParMatrix M;
        HypreParMatrix K;
        HypreParMatrix *T;    //T = M + dt K
        HypreParVector F; 
        CGSolver M_solver;    
        CGSolver T_solver;    
        HypreSmoother M_prec; 
        HypreSmoother T_prec; 

        //Extra
        mutable HypreParVector z;     //Auxiliar vector
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

        //Iteration parameters
        int iteration;
        double t;
        double dt;
        bool last;

        //Convergence parameters
        double actual_error;
        double total_error;
        int iterations_error;

        //Output parameters
        int dim;
        int serial_refinements;
        double h_min;
        HYPRE_Int size;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParGridFunction *x;
        HypreParVector *X;
        Array<int> ess_bdr; 
        FunctionCoefficient initial;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;

        //Extra
        ParaViewDataCollection *paraview_out;
};

extern double initial_conditions(const Vector &x, double t);

extern double T_f;     //Fusion temperature
extern double T_i;     //Initial temperature

extern double alpha_l; //Liquid thermal conduction
extern double alpha_s; //Solid thermal conduction
extern double lambda;   //s(t) = sqrt(4*lamda*(alpha_s+alpha_l)*t)
