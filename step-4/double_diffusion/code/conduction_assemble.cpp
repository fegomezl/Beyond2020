#include "header.h"

double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    TimeDependentOperator(2*fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    m_theta(NULL), m_phi(NULL),
    k_theta(NULL), k_phi(NULL),
    t_theta(NULL), t_phi(NULL),
    M_theta_solver(fespace.GetComm()), M_phi_solver(fespace.GetComm()),
    T_theta_solver(fespace.GetComm()), T_phi_solver(fespace.GetComm()),
    aux_theta(&fespace), aux_phi(&fespace), 
    aux_C(&fespace), aux_K(&fespace), aux_D(&fespace),
    aux_L(&fespace), aux_F(&fespace),
    psi(&fespace),
    coeff_r(r_f), zero(dim, zero_f), rot(dim, rot_f), 
    gradpsi(&psi), coeff_rV(rot, zero), dt_coeff_rV(0., coeff_rV),
    coeff_rCL(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r), dt_coeff_rK(0., coeff_rK),
    coeff_rD(coeff_r, coeff_r), dt_coeff_rD(0., coeff_rD),
    coeff_rCLV(coeff_r, zero), dt_coeff_rCLV(0., zero),
    rF_gradtheta(coeff_r, zero), dt_rF_gradtheta(0, zero)
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

    Array<int> ess_bdr_theta(attributes);
    ess_bdr_theta[0] = 1; ess_bdr_theta[1] = 1;
    ess_bdr_theta[2] = 0; ess_bdr_theta[3] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_list_theta);

    Array<int> ess_bdr_phi(attributes);
    ess_bdr_phi[0] = 1; ess_bdr_phi[1] = 1;
    ess_bdr_phi[2] = 0; ess_bdr_phi[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_list_phi);

    //Apply initial conditions
    ParGridFunction theta(&fespace);
    FunctionCoefficient initial_theta(initial_theta_f);
    theta.ProjectCoefficient(initial_theta);
    theta.ProjectBdrCoefficient(initial_theta, ess_bdr_theta);
    theta.GetTrueDofs(X.GetBlock(0));

    ParGridFunction phi(&fespace);
    FunctionCoefficient initial_phi(initial_phi_f);
    phi.ProjectCoefficient(initial_phi);
    phi.ProjectBdrCoefficient(initial_phi, ess_bdr_phi);
    phi.GetTrueDofs(X.GetBlock(1));

    //Configure M solvers
    M_theta_solver.iterative_mode = false;
    M_theta_solver.SetRelTol(config.reltol_conduction);
    M_theta_solver.SetAbsTol(config.abstol_conduction);
    M_theta_solver.SetMaxIter(config.iter_conduction);
    M_theta_solver.SetPrintLevel(0);
    M_theta_solver.SetPreconditioner(M_theta_prec);
    M_theta_prec.SetType(HypreSmoother::Jacobi);

    M_phi_solver.iterative_mode = false;
    M_phi_solver.SetRelTol(config.reltol_conduction);
    M_phi_solver.SetAbsTol(config.abstol_conduction);
    M_phi_solver.SetMaxIter(config.iter_conduction);
    M_phi_solver.SetPrintLevel(0);
    M_phi_solver.SetPreconditioner(M_phi_prec);
    M_phi_prec.SetType(HypreSmoother::Jacobi);

    //Configure T solvers
    T_theta_solver.iterative_mode = false;
    T_theta_solver.SetRelTol(config.reltol_conduction);
    T_theta_solver.SetAbsTol(config.abstol_conduction);
    T_theta_solver.SetMaxIter(config.iter_conduction);
    T_theta_solver.SetPrintLevel(0);
    T_theta_solver.SetPreconditioner(T_theta_prec);

    T_phi_solver.iterative_mode = false;
    T_phi_solver.SetRelTol(config.reltol_conduction);
    T_phi_solver.SetAbsTol(config.abstol_conduction);
    T_phi_solver.SetMaxIter(config.iter_conduction);
    T_phi_solver.SetPrintLevel(0);
    T_phi_solver.SetPreconditioner(T_phi_prec);

    SetParameters(X);
}

//Initial conditions
double initial_theta_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = Zmax/2;
    double sigma = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;
}

double initial_phi_f(const Vector &x){
    return 10;
}
