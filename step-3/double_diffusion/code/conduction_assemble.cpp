#include "header.h"

double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    m_theta(NULL), m_phi(NULL),
    k_theta(NULL), k_phi(NULL),
    t_theta(NULL), t_phi(NULL),
    M_theta(NULL), M_phi(NULL), M(NULL),
    K_theta(NULL), K_phi(NULL), K(NULL),
    T_theta(NULL), T_phi(NULL), T(NULL),
    M_blocks(2,2), K_blocks(2,2), T_blocks(2,2),
    M_solver(fespace.GetComm()), T_solver(fespace.GetComm()),
    SLU_M(NULL), SLU_T(NULL),
    aux_theta(&fespace), aux_phi(&fespace), 
    aux_C(&fespace), aux_K(&fespace), aux_D(&fespace),
    psi(&fespace),
    coeff_r(r_f), zero(dim, zero_f), rot(dim, rot_f), 
    gradpsi(&psi), coeff_rV(rot, zero), dt_coeff_rV(0., coeff_rV),
    coeff_rCL(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r), dt_coeff_rK(0., coeff_rK),
    coeff_rD(coeff_r, coeff_r), dt_coeff_rD(0., coeff_rD),
    coeff_rCLV(coeff_r, zero), dt_coeff_rCLV(0., zero)
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

    Array<int> ess_tdof_list_theta;
    Array<int> ess_bdr_theta(attributes);
    ess_bdr_theta[0] = 1; ess_bdr_theta[1] = 1;
    ess_bdr_theta[2] = 0; ess_bdr_theta[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_list_theta);

    Array<int> ess_tdof_list_phi;
    Array<int> ess_bdr_phi(attributes);
    ess_bdr_phi[0] = 1; ess_bdr_phi[1] = 1;
    ess_bdr_phi[2] = 0; ess_bdr_phi[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_list_phi);

    //Apply initial conditions
    ParGridFunction theta(&fespace);
    FunctionCoefficient initial_theta(initial_theta_f);
    theta.ProjectCoefficient(initial_theta);
    theta.ProjectBdrCoefficient(initial_theta, ess_bdr_theta);
    theta.ParallelAverage(X.GetBlock(0));

    ParGridFunction phi(&fespace);
    FunctionCoefficient initial_phi(initial_phi_f);
    phi.ProjectCoefficient(initial_phi);
    phi.ProjectBdrCoefficient(initial_phi, ess_bdr_phi);
    phi.ParallelAverage(X.GetBlock(1));

    //Configure M solver
    //M_solver.iterative_mode = false;
    M_solver.SetPrintStatistics(false);
    M_solver.SetSymmetricPattern(true);
    M_solver.SetColumnPermutation(superlu::PARMETIS);
    M_solver.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Configure T solver
    //T_solver.iterative_mode = false;
    T_solver.SetPrintStatistics(false);
    T_solver.SetSymmetricPattern(true);
    T_solver.SetColumnPermutation(superlu::PARMETIS);
    T_solver.SetIterativeRefine(superlu::SLU_DOUBLE);

    SetParameters(X);
}

//Initial conditions
double initial_theta_f(const Vector &x){
    return 0.001*(x(1) - Zmin)*(Zmax - x(1))*(x(0) + Rmax - 2*Rmin)*(Rmax - x(0));
}

double initial_phi_f(const Vector &x){
    return 0.002*(x(1) - Zmin)*(Zmax - x(1))*(x(0) + Rmax - 2*Rmin)*(Rmax - x(0));
}
