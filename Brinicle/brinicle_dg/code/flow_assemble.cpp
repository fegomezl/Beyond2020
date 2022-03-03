#include "header.h"

//Boundary values
double boundary_w(const Vector &x);
double boundary_psi(const Vector &x);

double boundary_psi_in(const Vector &x);
double boundary_psi_out(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_dg, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets):
    config(config),
    fespace(fespace), fespace_dg(fespace_dg), fespace_v(fespace_v),
    block_true_offsets(block_true_offsets),
    W(&fespace), Psi(&fespace),
    B_w(&fespace), B_psi(&fespace),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    M_e(NULL), D_e(NULL), C_e(NULL), Ct_e(NULL),
    psi(&fespace), w(&fespace),
    coeff_r(r_f), r_inv_hat(dim, r_inv_hat_f), rot(dim, rot_f), inv_R(inv_r), 
    grad(&fespace, &fespace_v)
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

    Array<int> ess_bdr_w(attributes);   ess_bdr_w = 0;
    ess_bdr_w[0] = 0; ess_bdr_w[1] = 0;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 0;
    ess_bdr_w[4] = 0; ess_bdr_w[5] = 0;
    fespace.GetEssentialTrueDofs(ess_bdr_w, ess_tdof_w);
  
    Array<int> ess_bdr_psi(attributes);   ess_bdr_psi = 0;
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    ess_bdr_psi[4] = 1; ess_bdr_psi[5] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_psi);

    Array<int> bdr_psi_in(attributes);   bdr_psi_in = 0;
    bdr_psi_in[4] = ess_bdr_psi[4];
    
    Array<int> bdr_psi_out(attributes);   bdr_psi_out = 0;
    bdr_psi_out[5] = ess_bdr_psi[5];
    
    Array<int> bdr_psi_closed_down(attributes);   bdr_psi_closed_down = 0;
    bdr_psi_closed_down[0] = ess_bdr_psi[0];
    bdr_psi_closed_down[2] = ess_bdr_psi[2];
    
    Array<int> bdr_psi_closed_up(attributes);   bdr_psi_closed_up = 0;
    bdr_psi_closed_up[1] = ess_bdr_psi[1];
    bdr_psi_closed_up[3] = ess_bdr_psi[3];

    //Setup boundary coefficients
    FunctionCoefficient w_coeff(boundary_w);
    FunctionCoefficient psi_coeff(boundary_psi);
    FunctionCoefficient psi_in(boundary_psi_in);
    FunctionCoefficient psi_out(boundary_psi_out);
    ConstantCoefficient closed_down(0.);
    ConstantCoefficient closed_up(Q);

    //Apply boundary conditions
    w.ProjectCoefficient(w_coeff);
    w.ProjectBdrCoefficient(closed_down, ess_bdr_w);
    w.ParallelProject(W);
 
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
    g.Assemble();
    g.ParallelAssemble(B_w);
    C_e->Mult(Psi, B_w, -1., 1.);
    EliminateBC(*M, *M_e, ess_tdof_w, W, B_w);

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
    double x_rel = x(0)/R_in;
    double y_rel = x(1)/Z_out;
    double in = 1., out = 1.;
    if (x(0) < R_in)
        in = pow(x_rel, 2)*(2-pow(x_rel, 2));
    if (x(1) < Z_out)
        out = pow(y_rel, 2)*(3-2*y_rel);
    return Q*in*out;
}

double boundary_psi_in(const Vector &x){
    double x_rel = x(0)/R_in;
    return Q*pow(x_rel, 2)*(2-pow(x_rel, 2));
}

double boundary_psi_out(const Vector &x){
    double y_rel = x(1)/Z_out;
    return Q*pow(y_rel, 2)*(3-2*y_rel);
}
