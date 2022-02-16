#include "header.h"

//Initial condition
double enthalphy_0(const Vector &x);

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Vector &X):
    config(config),
    TimeDependentOperator(fespace_L2.GetTrueVSize(), 0.),
    fespace_L2(fespace_L2), fespace_ND(fespace_ND),
    ess_bdr(attributes),
    m(NULL), k(NULL), b(NULL), 
    M(NULL), K(NULL), T(NULL),
    B(NULL), B_dt(NULL),
    Z(&fespace_L2),
    coeff_r(r_f), 
    coeff_zero(dim, zero_f),
    M_solver(MPI_COMM_WORLD), T_solver(MPI_COMM_WORLD)
{
    //Set boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0
    
    ess_bdr = 0;
    ess_bdr[0] = 0;  ess_bdr[1] = 0;
    ess_bdr[2] = 0;  ess_bdr[3] = 0;

    //Define initial condition
    ParGridFunction enthalphy(&fespace_L2);
    FunctionCoefficient coeff_enthalphy_0(enthalphy_0);
    enthalphy.ProjectCoefficient(coeff_enthalphy_0);
    enthalphy.GetTrueDofs(X);

    //Create mass matrix
    m = new ParBilinearForm(&fespace_L2);
    m->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();

    //Configure M solver
    M_prec.SetPrintLevel(0);
    M_prec.SetOperator(*M);

    M_solver.SetTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(*M);

    //Configure T solver
    T_prec.SetPrintLevel(0);

    T_solver.SetRelTol(config.reltol_conduction);
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction);
    T_solver.SetKDim(10);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);
}

double enthalphy_0(const Vector &x){
    return 1.;
}
