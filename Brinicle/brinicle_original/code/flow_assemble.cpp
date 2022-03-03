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
    Y(block_true_offsets), B(block_true_offsets),
    ess_bdr_psi(attributes), ess_bdr_w(attributes),
    bdr_psi_in(attributes), bdr_psi_out(attributes),
    bdr_psi_closed_down(attributes), bdr_psi_closed_up(attributes),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL), ct(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    psi(&fespace), w(&fespace), v(&fespace_v), rv(&fespace_v),
    theta(&fespace), theta_dr(&fespace), 
    phi(&fespace), phi_dr(&fespace), 
    phase(&fespace), eta(&fespace), psi_grad(&fespace_v), 
    coeff_r(r_f), inv_R(r_inv_f), r_inv_hat(dim, r_inv_hat_f),
    w_coeff(boundary_w), psi_coeff(boundary_psi), 
    psi_in(boundary_psi_in), psi_out(boundary_psi_out),
    closed_down(0.), closed_up(InflowFlux),
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
  
    ess_bdr_psi = 0;
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    ess_bdr_psi[4] = 1; ess_bdr_psi[5] = 1;

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
    w.ProjectCoefficient(w_coeff);
    w.ProjectBdrCoefficient(closed_down, ess_bdr_w);
 
    psi.ProjectCoefficient(psi_coeff);
    psi.ProjectBdrCoefficient(psi_in, bdr_psi_in);
    psi.ProjectBdrCoefficient(psi_out, bdr_psi_out);
    psi.ProjectBdrCoefficient(closed_down, bdr_psi_closed_down);
    psi.ProjectBdrCoefficient(closed_up, bdr_psi_closed_up);

    //Define the constant RHS
    g = new ParLinearForm(&fespace);
    g->Assemble();

    //Define constant bilinear forms of the system
    m = new ParBilinearForm (&fespace);
    m->AddDomainIntegrator(new MassIntegrator);
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, w, *g, Operator::DIAG_ONE);
    m->Finalize();
    M = m->ParallelAssemble();

    c = new ParMixedBilinearForm (&fespace, &fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator);
    c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c->Assemble();
    c->EliminateTrialDofs(ess_bdr_psi, psi, *g);
    c->EliminateTestDofs(ess_bdr_w);
    c->Finalize();
    C = c->ParallelAssemble();

    //Transfer to TrueDofs
    w.ParallelAssemble(Y.GetBlock(0));
    psi.ParallelAssemble(Y.GetBlock(1));
    g->ParallelAssemble(B.GetBlock(0));

    //Create gradient interpolator
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();

    //Set initial system
    SetParameters(X);
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
