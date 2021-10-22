#include "header.h"

//Initial conditions
double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    TimeDependentOperator(2*fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    m_theta(NULL), m_phi(NULL),
    k_theta(NULL), k_phi(NULL),
    M_theta(NULL), M_e_theta(NULL), M_0_theta(NULL), M_phi(NULL), M_e_phi(NULL), M_0_phi(NULL),
    K_0_theta(NULL),                                 K_0_phi(NULL),
    T_theta(NULL), T_e_theta(NULL),                  T_phi(NULL), T_e_phi(NULL),
    Z_theta(&fespace), Z_phi(&fespace),
    M_theta_solver(fespace.GetComm()), M_phi_solver(fespace.GetComm()),
    T_theta_solver(fespace.GetComm()), T_phi_solver(fespace.GetComm()),
    theta(&fespace), phi(&fespace), phase(&fespace), 
    aux_C(&fespace), aux_K(&fespace), aux_D(&fespace), aux_L(&fespace),
    coeff_r(r_f), zero(dim, zero_f), 
    coeff_rC(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r),
    coeff_rD(coeff_r, coeff_r),
    dHdT(zero, zero), dT_2(zero, zero)
{
    //Define essential boundary conditions
    //   
    //                  1
    //            /------------\
    //            |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0

    Array<int> ess_bdr_theta(attributes); ess_bdr_theta = 0;
    ess_bdr_theta  [0] = 1;   ess_bdr_theta  [1] = 1;
    ess_bdr_theta  [2] = 0;   ess_bdr_theta  [3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_theta);

    Array<int> ess_bdr_phi(attributes); ess_bdr_phi = 0;
    ess_bdr_phi  [0] = 0;     ess_bdr_phi  [1] = 0;
    ess_bdr_phi  [2] = 0;     ess_bdr_phi  [3] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_phi);

    //Apply initial conditions
    ParGridFunction theta(&fespace);
    FunctionCoefficient initial_theta(initial_theta_f);
    theta.ProjectCoefficient(initial_theta);
    theta.GetTrueDofs(X.GetBlock(0));

    ParGridFunction phi(&fespace);
    FunctionCoefficient initial_phi(initial_phi_f);
    phi.ProjectCoefficient(initial_phi);
    phi.GetTrueDofs(X.GetBlock(1));

    //Configure M solvers
    M_theta_solver.iterative_mode = false;
    M_theta_solver.SetTol(config.reltol_conduction);
    M_theta_solver.SetAbsTol(config.abstol_conduction);
    M_theta_solver.SetMaxIter(config.iter_conduction);
    M_theta_solver.SetPrintLevel(0);
    M_theta_prec.SetPrintLevel(0);
    M_theta_solver.SetPreconditioner(M_theta_prec);

    M_phi_solver.iterative_mode = false;
    M_phi_solver.SetTol(config.reltol_conduction);
    M_phi_solver.SetAbsTol(config.abstol_conduction);
    M_phi_solver.SetMaxIter(config.iter_conduction);
    M_phi_solver.SetPrintLevel(0);
    M_phi_prec.SetPrintLevel(0);
    M_phi_solver.SetPreconditioner(M_phi_prec);

    //Configure T solvers
    T_theta_solver.iterative_mode = false;
    T_theta_solver.SetTol(config.reltol_conduction);
    T_theta_solver.SetAbsTol(config.abstol_conduction);
    T_theta_solver.SetMaxIter(config.iter_conduction);
    T_theta_solver.SetPrintLevel(0);
    T_theta_prec.SetPrintLevel(0);
    T_theta_solver.SetPreconditioner(T_theta_prec);

    T_phi_solver.iterative_mode = false;
    T_phi_solver.SetTol(config.reltol_conduction);
    T_phi_solver.SetAbsTol(config.abstol_conduction);
    T_phi_solver.SetMaxIter(config.iter_conduction);
    T_phi_solver.SetPrintLevel(0);
    T_phi_prec.SetPrintLevel(0);
    T_phi_solver.SetPreconditioner(T_phi_prec);

    SetParameters(X);
}

//Initial conditions
double initial_theta_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = (Zmax + Zmin)/2;
    double Rad = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(Rad, 2))
        return -5.;
    else
        return -1.8;
}

double initial_phi_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    if (x(0) < mid_x)
        return 10.;
    else
        return 3.5;
}
