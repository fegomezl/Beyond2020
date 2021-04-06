#pragma once

#include <iostream>
#include <fstream>
#include <string>
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
        Conduction_Operator(ParFiniteElementSpace *&fespace, double t_init,
                            double alpha, double kappa, const Vector &X);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt, 
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        void SetParameters(const Vector &X);                        //Update the bilinear forms

        virtual ~Conduction_Operator();
    protected:
        //Operator parameters
        double t_init;
        double current_dt;
        double alpha;
        double kappa;

        //Mesh objects
        ParFiniteElementSpace *fespace;
        Array<int> ess_tdof_list; 

        //System objects
        ParBilinearForm *m;  //Mass operator
        ParBilinearForm *k;  //Difussion operator

        //Solver objects
        HypreParMatrix M;
        HypreParMatrix K;
        HypreParMatrix *T;    //T = M + dt K
        CGSolver M_solver;    
        CGSolver T_solver;    
        HypreSmoother M_prec; 
        HypreSmoother T_prec; 

        //Extra
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
        bool delete_fec;   //Need to delete fec at the end
};

double initial_conditions(const Vector &X);

extern double int_rad;
extern double out_rad;
