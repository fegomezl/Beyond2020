#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "mfem.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

using namespace std;
using namespace mfem;

struct Config{
    Config(bool master, int nproc);

    bool master;
    int nproc;
    int order;
    int refinements;
    double dt_init;
    double t_final;
    int vis_steps_max;
    int ode_solver_type;
    double reltol;
    double abstol;
    double print;
};

class Conduction_Operator : public TimeDependentOperator{
    public:
        Conduction_Operator(ParFiniteElementSpace &fespace, Array<int> ess_bdr);

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Solver for explicit methods
        virtual void ImplicitSolve(const double dt,
                                   const Vector &X, Vector &dX_dt); //Solver for implicit methods
        virtual int SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt);
   	    virtual int SUNImplicitSolve(const Vector &B, Vector &X, double tol);

        virtual ~Conduction_Operator();

    protected:
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

        FunctionCoefficient r;
        ProductCoefficient r_alpha;
        ProductCoefficient dt_r_alpha;
};

class Artic_sea{
    public:
        Artic_sea(Config config);
        void run(const char *mesh_file, double &H_size, double &Total_time, double &Mean_error);
        ~Artic_sea();
    private:
        void make_grid(const char *mesh_file);
        void assemble_system();
        void time_step();

        //Global parameters
        Config config;

        int iteration;
        double t;
        double dt;
        bool last;
        int vis_iteration;
        int vis_steps;
        int vis_impressions;
        double h_size;
        double actual_error;
  	    double total_error;

        int dim;
        int serial_refinements;
        HYPRE_Int size;

        ParMesh *pmesh;
        FiniteElementCollection *fec;
        ParFiniteElementSpace *fespace;

        //System objects
        ParGridFunction *x;
        HypreParVector *X;
        Array<int> ess_bdr;
        FunctionCoefficient initial_f;
        Conduction_Operator *oper;

        //Solver objects
        ODESolver *ode_solver;
        CVODESolver *cvode;
        ARKStepSolver *arkode;
};

extern double Rmin, Rmax, Zmin, Zmax;
extern int  Mterms;
extern int Nterms;
extern std::vector<double> Coeficients;

extern double alpha;

extern void Calc_Coe(double a, double b, std::vector<double> & Coeficients);
extern double initial(const Vector &x, double t);
extern double initial_condition(double r, double z, int m, int n, double a, double b);
extern double integrand(double r,double z,int m,int n, double a, double b);
extern double integrate(int m, int n,double a,double b);
extern double Aux( double r, double z, double t);
