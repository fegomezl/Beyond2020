#include "header.h"

//Boundary values
double boundary_w(const Vector &x);
double boundary_psi(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const Vector &Theta):
    config(config),
    fespace(fespace),
    block_true_offsets(3),
    ess_bdr_psi(attributes), ess_bdr_w(attributes),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL), ct(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    psi(&fespace), w(&fespace), v(&fespace_v), rv(&fespace_v),
    theta(&fespace), theta_eta(&fespace), 
    psi_grad(&fespace_v), theta_dr(&fespace),
    coeff_r(r_f), inv_R(inv_r), r_inv_hat(dim, r_inv_hat_f),
    w_coeff(boundary_w), psi_coeff(boundary_psi), 
    grad(&fespace, &fespace_v),
    rot(dim, rot_f), Psi_grad(&psi_grad),
    rot_Psi_grad(rot, Psi_grad)
{
    //Create the block offsets
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace.TrueVSize();
    block_true_offsets[2] = fespace.TrueVSize();
    block_true_offsets.PartialSum();
  
    //Initialize the corresponding vectors
    X.Update(block_true_offsets);
    B.Update(block_true_offsets);
  
    //Define essential boundary conditions
    //
    //                  1
    //            /------------\
    // (w,psi=0)  |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0
    //
    ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
  
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;

    //Apply boundary conditions
    w.ProjectCoefficient(w_coeff);
    psi.ProjectCoefficient(psi_coeff);

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
    w.ParallelAssemble(X.GetBlock(0));
    psi.ParallelAssemble(X.GetBlock(1));
    g->ParallelAssemble(B.GetBlock(0));

    //Create gradient interpolator
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();

    //Set initial system
    SetParameters(Theta);
}

//Boundary values
double boundary_w(const Vector &x){
    return 0.;
}

double boundary_psi(const Vector &x){
    return 0.;
}
