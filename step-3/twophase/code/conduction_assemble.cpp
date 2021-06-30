#include "header.h"

double initial_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int attributes, Vector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    m(NULL),
    k(NULL),
    t(NULL),
    aux(&fespace),
    aux_C(&fespace),
    aux_K(&fespace),
    r(rf),
    coeff_C(&aux_C), coeff_rC(r, coeff_C),
    coeff_K(&aux_K), coeff_rK(r, coeff_K), dt_coeff_rK(1., coeff_rK),
    coeff_L(&aux), coeff_rL(r, coeff_L),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm())
{
    const double rel_tol = 1e-8;

    //Set boundary conditions
    Array<int> ess_bdr(attributes);
    ess_bdr = 0;  ess_bdr[0] = 1;  ess_bdr[1] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    FunctionCoefficient initial(initial_f);

    //Define solution x and apply initial conditions
    ParGridFunction x(&fespace);
    x.ProjectCoefficient(initial);
    x.ProjectBdrCoefficient(initial, ess_bdr);
    x.GetTrueDofs(X);

    //Configure M solver
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_prec.SetType(HypreSmoother::Jacobi);

    //Configure T solver
    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(X);
}

double initial_f(const Vector &x){
    double mid = 0.6*Zmax;
    if (x(1) <= mid)
        return -10*(1 - x(1)/mid);
    else
        return 10*(x(1) - mid)/(Zmax - mid);
}
