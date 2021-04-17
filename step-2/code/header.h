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
    double t_final;
    int vis_steps;
    int ode_solver_type;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(ParFiniteElementSpace *&fespace, const Vector &X, double b_size);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt, 
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        void SetParameters(const Vector &X, Array<int> ess_bdr);                        //Update the bilinear forms

        virtual ~Conduction_Operator();
    protected:
        //Operator parameters
        double current_dt;

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

extern double d_bdr(const Vector &x, double t);

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
        int serial_refinements;
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
        FunctionCoefficient boundary = FunctionCoefficient(d_bdr);

        //Solver objects
        ODESolver *ode_solver;

        //Extra
        ParaViewDataCollection *paraview_out;
};

extern double T_f;     //Fusion temperature
extern double T_i;     //Initial temperature

extern double alpha_l; //Liquid thermal conduction
extern double alpha_s; //Solid thermal conduction
extern double lambda;   //s(t) = sqrt(4*lamda*(alpha_s+alpha_l)*t)