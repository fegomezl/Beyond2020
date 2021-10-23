#include "header.h"

//Initial conditions
double initial_theta_f(const Vector &x);
double initial_phi_f(const Vector &x);

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    TimeDependentOperator(2*fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    M(NULL), T(NULL),
    M_0(NULL), K_0(NULL), T_e(NULL),
    Z(block_true_offsets),
    dZ(block_true_offsets), T_d(block_true_offsets),
    M_solver(NULL), T_solver(NULL),
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

    Array<int> ess_tdof_theta;
    Array<int> ess_bdr_theta(attributes); ess_bdr_theta = 0;
    ess_bdr_theta  [0] = 1;   ess_bdr_theta  [1] = 1;
    ess_bdr_theta  [2] = 0;   ess_bdr_theta  [3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_theta, ess_tdof_theta);

    Array<int> ess_tdof_phi;
    Array<int> ess_bdr_phi(attributes); ess_bdr_phi = 0;
    ess_bdr_phi  [0] = 0;     ess_bdr_phi  [1] = 0;
    ess_bdr_phi  [2] = 0;     ess_bdr_phi  [3] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_phi, ess_tdof_phi);

    ess_tdof_list = ess_tdof_theta;
    ess_tdof_list.Append(ess_tdof_phi);
    for (int ii = ess_tdof_theta.Size(); ii < ess_tdof_list.Size(); ii++)
        ess_tdof_list[ii] += block_true_offsets[1];

    //Apply initial conditions
    ParGridFunction theta(&fespace);
    FunctionCoefficient initial_theta(initial_theta_f);
    theta.ProjectCoefficient(initial_theta);
    theta.GetTrueDofs(X.GetBlock(0));

    ParGridFunction phi(&fespace);
    FunctionCoefficient initial_phi(initial_phi_f);
    phi.ProjectCoefficient(initial_phi);
    phi.GetTrueDofs(X.GetBlock(1));

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
