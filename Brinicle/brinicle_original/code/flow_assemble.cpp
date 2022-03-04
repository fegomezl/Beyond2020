#include "header.h"

//Boundary values
double boundary_w(const Vector &x);
double boundary_psi(const Vector &x);

double boundary_psi_in(const Vector &x);
double boundary_psi_out(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, BlockVector &X):
    config(config),
    fespace(fespace),
    block_true_offsets(block_true_offsets),
    ess_bdr_psi(attributes), ess_bdr_w(attributes),
    bdr_psi_in(attributes), bdr_psi_out(attributes),
    bdr_psi_closed_down(attributes), bdr_psi_closed_up(attributes),
    W(&fespace), Psi(&fespace),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    F(&fespace), G(&fespace), B(block_true_offsets),
    M_e(NULL), D_e(NULL), C_e(NULL), Ct_e(NULL),
    v(&fespace_v), rv(&fespace_v),
    theta(&fespace), theta_dr(&fespace), 
    phi(&fespace), phi_dr(&fespace), 
    phase(&fespace), eta(&fespace), 
    psi(&fespace),
    psi_grad(&fespace_v), 
    coeff_r(r_f), inv_R(r_inv_f), r_inv_hat(dim, r_inv_hat_f),
    grad(&fespace, &fespace_v),
    rot(dim, rot_f), Psi_grad(&psi_grad),
    rot_Psi_grad(rot, Psi_grad)
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
    //

    ess_bdr_w = 0;
    ess_bdr_w[0] = 0; ess_bdr_w[1] = 0;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 0;
    ess_bdr_w[4] = 0; ess_bdr_w[5] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_w, ess_tdof_w);
  
    ess_bdr_psi = 0;
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    ess_bdr_psi[4] = 1; ess_bdr_psi[5] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_psi);

    bdr_psi_in = 0;
    bdr_psi_in[4] = ess_bdr_psi[4];
    
    bdr_psi_out = 0;
    bdr_psi_out[5] = ess_bdr_psi[5];
    
    bdr_psi_closed_down = 0;
    bdr_psi_closed_down[0] = ess_bdr_psi[0];
    bdr_psi_closed_down[2] = ess_bdr_psi[2];
    
    bdr_psi_closed_up = 0;
    bdr_psi_closed_up[1] = ess_bdr_psi[1];
    bdr_psi_closed_up[3] = ess_bdr_psi[3];

    //Apply boundary conditions
    FunctionCoefficient w_coeff(boundary_w);
    ParGridFunction w(&fespace);
    w.ProjectCoefficient(w_coeff);
    w.ProjectBdrCoefficient(w_coeff, ess_bdr_w);
    w.ParallelProject(W);
 
    FunctionCoefficient psi_coeff(boundary_psi);
    FunctionCoefficient psi_in(boundary_psi_in);
    FunctionCoefficient psi_out(boundary_psi_out);
    ConstantCoefficient closed_down(0.);
    ConstantCoefficient closed_up(InflowFlux);
    ParGridFunction psi(&fespace);
    psi.ProjectCoefficient(psi_coeff);
    psi.ProjectBdrCoefficient(psi_in, bdr_psi_in);
    psi.ProjectBdrCoefficient(psi_out, bdr_psi_out);
    psi.ProjectBdrCoefficient(closed_down, bdr_psi_closed_down);
    psi.ProjectBdrCoefficient(closed_up, bdr_psi_closed_up);
    psi.ParallelProject(Psi);

    //Define constant bilinear forms of the system
    ParBilinearForm m(&fespace);
    m.AddDomainIntegrator(new MassIntegrator);
    m.Assemble();
    m.Finalize();
    M = m.ParallelAssemble();
    M_e = M->EliminateRowsCols(ess_tdof_w);

    ParMixedBilinearForm c(&fespace, &fespace);
    c.AddDomainIntegrator(new MixedGradGradIntegrator);
    c.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c.Assemble();
    c.Finalize();
    C = c.ParallelAssemble();
    C->EliminateRows(ess_tdof_w);
    C_e = C->EliminateCols(ess_tdof_psi);

    Ct = C->Transpose();
    Ct->EliminateRows(ess_tdof_psi);
    Ct_e = Ct->EliminateCols(ess_tdof_w);

    //Define the constant RHS
    ParLinearForm g(&fespace);
    ConstantCoefficient Zero(0.);
    g.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    g.Assemble();
    g.ParallelAssemble(B.GetBlock(0));
    C_e->Mult(Psi, G, -1., 1.);
    EliminateBC(*M, *M_e, ess_tdof_w, W, G);

    //Create gradient interpolator
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();
}

//Boundary values
double boundary_w(const Vector &x){
    return 0.;
}

double boundary_psi(const Vector &x){
    double x_rel = x(0)/RIn;
    double y_rel = x(1)/ZOut;
    double in = 1., out = 1.;
    if (x(0) < RIn)
        in = pow(x_rel, 2)*(2-pow(x_rel, 2));
    if (x(1) < ZOut)
        out = pow(y_rel, 2)*(3-2*y_rel);
    return InflowFlux*in*out;
}

double boundary_psi_in(const Vector &x){
    double x_rel = x(0)/RIn;
    return InflowFlux*pow(x_rel, 2)*(2-pow(x_rel, 2));
}

double boundary_psi_out(const Vector &x){
    double y_rel = x(1)/ZOut;
    return InflowFlux*pow(y_rel, 2)*(3-2*y_rel);
}
