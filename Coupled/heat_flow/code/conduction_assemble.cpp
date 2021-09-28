#include "header.h"

//Initial temperature field
double initial_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Vector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    m(NULL), k(NULL), 
    M(NULL), M_e(NULL), M_0(NULL),
    K_0(NULL), 
    T(NULL), T_e(NULL),
    Z(&fespace),
    aux(&fespace), aux_C(&fespace), aux_K(&fespace), aux_L(&fespace),
    rv(&fespace_v),
    coeff_r(r_f), zero(dim, zero_f),
    coeff_rC(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r), 
    coeff_rCV(coeff_r, zero),
    dHdT(zero, zero), dT_2(zero, zero),
    M_solver(fespace.GetComm()), T_solver(fespace.GetComm())
{
    //Set boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //    (0)    2|            |3
    //            |            |
    //            \------------/
    //                  0

    Array<int> ess_bdr(attributes);
    ess_bdr[0] = 1;  ess_bdr[1] = 1;
    ess_bdr[2] = 0;  ess_bdr[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Define initial condition
    FunctionCoefficient initial(initial_f);    
    aux.ProjectCoefficient(initial);
    aux.ProjectBdrCoefficient(initial, ess_bdr);
    aux.GetTrueDofs(X);

    //Configure M solver
    M_solver.SetTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_prec.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);

    //Configure T solver
    T_solver.SetTol(config.reltol_conduction);
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction);
    T_solver.SetPrintLevel(0);
    T_prec.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);
}

double initial_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = (Zmax + Zmin)/2;
    double Rad = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(Rad, 2))
        return -10;
    else
        return 10;
}
