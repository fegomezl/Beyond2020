#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Config{
    //Passing parameters
    bool master;
    int order;
    int serial_refinements;
    int refinements;
    int ode_solver_type;
    double t_init;
    double t_final;
    double dt_init;
    double alpha;
    double kappa;
    int vis_steps;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(ParFiniteElementSpace &f, 
                           double alpha, double kappa, const Vector &u);
        virtual void Mult(const Vector &u, Vector &du_dt) const; 
        //Solve the Backward-Euler equation k = f(u + dt*k, t) for k
        virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);
        //Update the diffusion BillinearForm K using the trueDOF vector X
        void SetParameters(const Vector &u);
        virtual ~Conduction_Operator();

    protected:
        ParFiniteElementSpace &fespace;
        Array<int> ess_tdof_list; 

        ParBilinearForm *M;
        ParBilinearForm *K;

        HypreParMatrix Mmat;
        HypreParMatrix Kmat;
        HypreParMatrix *T;    // T = M + dt K
        double current_dt;

        CGSolver M_Solver;    //Krylov solver for M^-1
        HypreSmoother M_prec; //Preconditioner por M

        CGSolver T_Solver;    //Implicit solver for T = M + dt K
        HypreSmoother T_prec; //Preconditioner por implicit solver

        double alpha, kappa;  //Nonlinear parameters

        mutable Vector z;     //Auxiliar vector
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
        int vis_steps;
        bool last;

        //Output parameters
        int dim;
        double h_min;
        HYPRE_Int size;

        //Mesh objects
        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParGridFunction *x;
        Vector X;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;

        //Extra
        ParaViewDataCollection *paraview_out;
        bool delete_fec;
};

double initial_conditions(const Vector &x);

extern double height;
extern double int_rad;
extern double out_rad;
