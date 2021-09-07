#include "header.h"

double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

double robin_h_theta_f(const Vector &x);
double robin_h_phi_f(const Vector &x);
double robin_ref_theta_f(const Vector &x);
double robin_ref_phi_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    TimeDependentOperator(2*fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    robin_bdr_theta(attributes), robin_bdr_phi(attributes),
    m_theta(NULL), m_phi(NULL),
    k_theta(NULL), k_phi(NULL),
    t_theta(NULL), t_phi(NULL),
    M_theta_solver(fespace.GetComm()), M_phi_solver(fespace.GetComm()),
    T_theta_solver(fespace.GetComm()), T_phi_solver(fespace.GetComm()),
    aux_theta(&fespace), aux_phi(&fespace),
    aux_C(&fespace), aux_K(&fespace), aux_D(&fespace),
    v(&fespace_v),
    coeff_r(r_f), zero(dim, zero_f),
    coeff_rV(&v), dt_coeff_rV(0., coeff_rV),
    coeff_rCL(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r), dt_coeff_rK(0., coeff_rK),
    coeff_rD(coeff_r, coeff_r), dt_coeff_rD(0., coeff_rD),
    coeff_rCLV(coeff_r, zero),  dt_coeff_rCLV(0., zero),
    robin_h_theta(robin_h_theta_f),                  robin_h_phi(robin_h_phi_f),
    robin_ref_theta(robin_ref_theta_f),              robin_ref_phi(robin_ref_phi_f),
    r_robin_h_theta(coeff_r, robin_h_theta),         r_robin_h_phi(coeff_r, robin_h_phi),
    r_robin_h_ref_theta(r_robin_h_theta,             robin_ref_theta), r_robin_h_ref_phi(r_robin_h_phi, robin_ref_phi),
    dt_r_robin_h_theta(0., r_robin_h_theta),         dt_r_robin_h_phi(0., r_robin_h_phi),
    dt_r_robin_h_ref_theta(0., r_robin_h_ref_theta), dt_r_robin_h_ref_phi(0., r_robin_h_ref_phi)
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

    Array<int> ess_bdr_theta(attributes); ess_bdr_theta = 0;
    ess_bdr_theta[0] = 0;     ess_bdr_theta[1] = 0;     ess_bdr_theta[3] = 0;
    robin_bdr_theta[0] = 0;   robin_bdr_theta[1] = 1;   robin_bdr_theta[3] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_list_theta);

    Array<int> ess_bdr_phi(attributes); ess_bdr_phi = 0;
    ess_bdr_phi[0] = 0;       ess_bdr_phi[1] = 0;       ess_bdr_phi[3] = 0;
    robin_bdr_phi[0] = 0;     robin_bdr_phi[1] = 1;     robin_bdr_phi[3] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_list_phi);

    //Check that the internal boundaries is always zero.
    ess_bdr_theta[2] = 0;     ess_bdr_phi[2] = 0;
    robin_bdr_theta[2] = 0;   robin_bdr_phi[2] = 0;

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
}

//Initial conditions

double initial_theta_f(const Vector &x){
    return theta_in;
}

double initial_phi_f(const Vector &x){
    return phi_in;
}

//Robin boundary conditions of the form kdu/dn = h(u-u_ref)

double robin_h_theta_f(const Vector &x){
    if (x(0) < r0)
        return vel*c_l;
    else
        return 0;
}

double robin_h_phi_f(const Vector &x){
    if (x(0) < r0)
        return vel;
    else
        return 0;
}

double robin_ref_theta_f(const Vector &x){
    return theta_out;
}

double robin_ref_phi_f(const Vector &x){
    return phi_out;
}
