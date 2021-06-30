#include "header.h"

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr):
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

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

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
