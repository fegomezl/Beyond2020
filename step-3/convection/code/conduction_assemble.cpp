#include "header.h"

double initial_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Vector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    m(NULL), k(NULL), t(NULL),
    aux(&fespace), aux_C(&fespace), aux_K(&fespace),
    psi(&fespace),
    r(r_f), zero(dim, zero_f),
    coeff_rCL(r, r),
    coeff_rK(r, r), dt_coeff_rK(0., coeff_rK),
    rot(dim, rot_f), gradpsi(&psi), rV(rot, zero),
    coeff_rCLV(r, zero), dt_coeff_rCLV(0., zero),
    M_solver(fespace.GetComm()), T_solver(fespace.GetComm())
{
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
    M_solver.SetRelTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_prec.SetType(HypreSmoother::Jacobi);

    //Configure T solver
    T_solver.iterative_mode = false;
    T_solver.SetRelTol(config.reltol_conduction);
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction);
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
