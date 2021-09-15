#include "header.h"

//Initial temperature field
double initial_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Vector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    m(NULL), k(NULL), t(NULL),
    aux(&fespace), aux_C(&fespace), aux_K(&fespace),
    psi(&fespace), v(&fespace_v),
    r(r_f), zero(dim, zero_f),
    coeff_rCL(r, r),
    coeff_rK(r, r), dt_coeff_rK(0., coeff_rK),
    rV(&v),
    coeff_rCLV(r, zero), dt_coeff_rCLV(0., zero),
    M_solver(fespace.GetComm()), T_solver(fespace.GetComm())
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
    Array<int> ess_bdr(attributes);
    ess_bdr[0] = 0;  ess_bdr[1] = 0;
    ess_bdr[2] = 0;  ess_bdr[3] = 0;
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
}

double initial_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = (Zmax - Zmin)/2;
    double Rad = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(Rad, 2))
        return 10;
    else
        return 20;
}
