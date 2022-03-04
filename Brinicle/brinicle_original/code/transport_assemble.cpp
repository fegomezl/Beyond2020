#include "header.h"

//Initial conditions
double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    config(config),
    TimeDependentOperator(2*fespace.GetTrueVSize(), config.t_init),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    m_theta(NULL), m_phi(NULL),
    k_theta(NULL), k_phi(NULL),
    M_theta(NULL), M_e_theta(NULL), M_0_theta(NULL), M_phi(NULL), M_e_phi(NULL), M_0_phi(NULL),
    K_0_theta(NULL),                                 K_0_phi(NULL),
    T_theta(NULL), T_e_theta(NULL),                  T_phi(NULL), T_e_phi(NULL),
    F_theta(&fespace), F_phi(&fespace),
    dt_F_theta(&fespace), dt_F_phi(&fespace),
    Z_theta(&fespace), Z_phi(&fespace),
    M_theta_solver(fespace.GetComm()), M_phi_solver(fespace.GetComm()),
    T_theta_solver(fespace.GetComm()), T_phi_solver(fespace.GetComm()),
    theta(&fespace), phi(&fespace), phase(&fespace), rv(&fespace_v), 
    aux_C(&fespace), aux_K(&fespace), aux_D(&fespace), aux_L(&fespace),
    coeff_r(r_f), zero(dim, zero_f), 
    coeff_rC(coeff_r, coeff_r),
    coeff_rK(coeff_r, coeff_r),
    coeff_rD(coeff_r, coeff_r),
    coeff_rV(&rv), coeff_rCV(coeff_r, zero), 
    dHdT(zero, zero), dT_2(zero, zero)
{
    //Define essential boundary conditions
    //   
    //              4               1
    //            /---|---------------------------\
    //            |                               |
    //            |                               |
    //            |                               | 3
    //            |                               |
    //          2 |                               |
    //            |                               -
    //            |                               |
    //            |                               | 5
    //            |                               |
    //            \-------------------------------/
    //                            0

    Array<int> ess_bdr_theta(attributes); ess_bdr_theta = 0;
    ess_bdr_theta  [0] = 0;   ess_bdr_theta  [1] = 0;   ess_bdr_theta  [3] = 0;
    ess_bdr_theta  [4] = 1;   ess_bdr_theta  [5] = 0;   
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_theta);

    Array<int> ess_bdr_phi(attributes); ess_bdr_phi = 0;
    ess_bdr_phi  [0] = 0;     ess_bdr_phi  [1] = 0;     ess_bdr_phi  [3] = 0; 
    ess_bdr_phi  [4] = 1;     ess_bdr_phi  [5] = 0;      
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_phi);

    //Check that the internal boundaries is always zero.
    ess_bdr_theta  [2] = 0;   ess_bdr_phi  [2] = 0; 

    //Apply initial conditions
    if (!config.restart){
        ParGridFunction theta(&fespace);
        FunctionCoefficient initial_theta(initial_theta_f);
        ConstantCoefficient theta_nu(NucleationTemperature);
        theta.ProjectCoefficient(initial_theta);
        theta.ProjectBdrCoefficient(theta_nu, ess_bdr_theta);
        theta.GetTrueDofs(X.GetBlock(0));
        
        ParGridFunction phi(&fespace);
        FunctionCoefficient initial_phi(initial_phi_f);
        ConstantCoefficient phi_nu(NucleationSalinity);
        phi.ProjectCoefficient(initial_phi);
        phi.ProjectBdrCoefficient(phi_nu, ess_bdr_phi);
        phi.GetTrueDofs(X.GetBlock(1));
    }

    //Set RHS
    ParLinearForm f_theta(&fespace);
    f_theta.Assemble();
    f_theta.ParallelAssemble(F_theta);
    F_theta = 0.;

    ParLinearForm f_phi(&fespace);
    f_phi.Assemble();
    f_phi.ParallelAssemble(F_phi);
    F_phi = 0.;

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
}

//Initial conditions

double initial_theta_f(const Vector &x){
    if ((x(0) > RIn && x(1) > ZMax - NucleationHeight) || pow(x(0)-(RIn+ NucleationLength), 2) + pow(x(1)-(ZMax - NucleationHeight), 2) < pow(NucleationLength, 2))
        return -10;
    else
        return InitialTemperature;
}

double initial_phi_f(const Vector &x){
  // if (x(0) < R_in && x(1) > Zmax - n_h/4)
  // return ((phi_out-phi_in)/n_h)*(x(1)-Zmax+n_h*1.2)+phi_in;
      
    if ((x(0) > RIn && x(1) > ZMax - NucleationHeight) || pow(x(0)-(RIn+ NucleationLength), 2) + pow(x(1)-(ZMax - NucleationHeight), 2) < pow(NucleationLength, 2))
        return InitialSalinity;
    else
        return InitialSalinity;
}
