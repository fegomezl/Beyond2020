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
};

class Transport_Operator : public TimeDependentOperator{        //Evolution of the conserved quantities (entalphy and salinity)
    public:
        Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X);

        void SetParameters(const BlockVector &X, const Vector &r_Velocity);    //Update parameters from previous step

        virtual void Mult(const Vector &X, Vector &dX_dt) const;    //Standard solver
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);  //Sundials setup
	    virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);   //Sundials solver

        virtual ~Transport_Operator();
    protected:
        //Global parameters
        Config config;
 
        ParFiniteElementSpace &fespace_H1;
        ParFiniteElementSpace &fespace_ND;
        Array<int> block_offsets_H1;
        Array<int> ess_tdof_0, ess_tdof_1;

        //System objects
        ParBilinearForm *m0, *m1;  //Mass operator
        ParBilinearForm *k0, *k1;  //Difussion and convection operator
        ParLinearForm   *b0, *b1;  //Sources and boundary terms

        HypreParMatrix *M0,   *M1;  
        HypreParMatrix *M0_o, *M1_o;
        HypreParMatrix *M0_e, *M1_e;
        HypreParMatrix *K0,   *K1;
        HypreParMatrix *T0,   *T1;
        HypreParMatrix *T0_e, *T1_e;

        HypreParVector *B0,    *B1;
        HypreParVector *B0_dt, *B1_dt;

        mutable HypreParVector Z0, Z1;   

        HypreBoomerAMG M0_prec, M1_prec;
        HypreBoomerAMG T0_prec, T1_prec;
        HyprePCG M0_solver, M1_solver;
        HyprePCG T0_solver, T1_solver;

        FunctionCoefficient coeff_r;
        VectorFunctionCoefficient coeff_zero;
};

class Flow_Operator{
    public:
        Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1);

        void SetParameters(const BlockVector &X);

        void Solve(BlockVector &Y, Vector &Velocity, Vector &r_Velocity);

        ~Flow_Operator();
    protected:
        //Global parameters
        Config config;

        //Mesh objects
        ParFiniteElementSpace &fespace_H1;
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
        FiniteElementCollection *fec_ND;

        ParFiniteElementSpace *fespace_H1;
        ParFiniteElementSpace *fespace_ND;

        Array<int> block_offsets_H1;

        int dim;
        double h_min;
        int serial_refinements;

        HYPRE_Int size_H1;
        HYPRE_Int size_ND;

        //System objects
        ParGridFunction *enthalphy;
        ParGridFunction *salinity;
        ParGridFunction *vorticity;
        ParGridFunction *stream;
        ParGridFunction *velocity;
        ParGridFunction *r_velocity;

        BlockVector X;          //Enthalphy and Salinity
        BlockVector Y;          //Vorticity and Stream
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
extern double R, Z;         //Dimensions of the domain
extern double R_in, Z_out;  //Dimensions of the flow entrance and exit
extern double LenghtScale;  //Scaling for the lenght dimensions
extern double TimeScale;    //Scaling for the time dimensions
extern double JumpScale;    //Size of the indetermination window in heaviside functions (inverse)
extern double Epsilon;      //Size of the indetermination window in heaviside functions

extern double InflowVelocity;           //Velocity of the inflow
extern double InflowFlux;               //Flux of the given inflow velocity field
extern double InitialTemperature;       //Initial temperature of the domain
extern double InflowTemperature;        //Temperature of the inflow
extern double InitialSalinity;          //Initial salinity of the domain
extern double InflowSalinity;           //Salinity of the inflow
extern double NucleationLength;         //Lenght of the nucleation point
extern double NucleationHeight;         //Height of the nucleation point
extern double NucleationTemperature;    //Temperature of the nucleation point
extern double NucleationSalinity;       //Salinity of the nucleation point

//Rotational functions
extern double r_f(const Vector &x);                     //Function for r
extern double r_inv_f(const Vector &x);                 //Function for 1/r
extern void zero_f(const Vector &x, Vector &f);         //Function for 0 (vector)
extern void r_inv_hat_f(const Vector &x, Vector &f);    //Function for (1/r)*r^ (r^ unitary vector)
extern void rot_f(const Vector &x, DenseMatrix &f);     //Function for ( 0   1 )
                                                        //             (-1   0 )

//Physical properties (in H,S)
extern double FusionPoint(const double S);                          //Fusion temperature at a given salinity
extern double Phase(const double H);                                //Phase indicator (1 for liquid and 0 for solid)
extern double HeatDiffusivity(const double H);                      //Diffusion coefficient for the heat equation
extern double SaltDiffusivity(const double H);                      //Diffusion coefficient for the mass equation
extern double FusionDiffusivity(const double H, const double S);    //Diffusion coefficient for the latent term (-kGrad(T_f))
extern double InversePermeability(const double H);                  //Inverse of the brinkman penalization permeability
extern double ExpansivityEnthalphy(const double H, const double S); //Expansivity coefficient for the enthalphy gradient
extern double ExpansivitySalinity(const double H, const double S);  //Expansivity coefficient for the salinity gradient
extern double Buoyancy(const double H, const double S);             //Buoyancy coefficient

//Relationship between variables
extern double TStoHS(const double T, const double S);  //Map from temperature-salinity to enthaphy-salinity
extern double HStoTS(const double H, const double S);  //Map from enthalphy-salinity to temperature-salinity
