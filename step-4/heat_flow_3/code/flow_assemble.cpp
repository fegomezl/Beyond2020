#include "header.h"

//Unitary r vector
void r_hat_f(const Vector &x, Vector &y);

//Inverse of r
double r_inv_f(const Vector &x);

//Boundary values for w
double boundary_w(const Vector &x);

void boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
double boundary_psi(const Vector &x);

void boundary_gradpsi(const Vector &x, Vector &f);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, const HypreParVector *Theta):
    config(config),
    fespace(fespace),
    block_true_offsets(3),
    ess_bdr_psi(attributes), ess_bdr_w(attributes),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL),
    M(NULL), D(NULL), C(NULL),
    psi(NULL), w(NULL), v(NULL),
    w_aux(NULL), psi_aux(NULL), v_aux(NULL),
    theta_aux(NULL), theta_grad_aux(NULL), theta_eta(NULL), theta_rho(NULL),
    r(r_f), r_inv(r_inv_f),
    r_hat(dim, r_hat_f), r_inv_hat(r_inv, r_hat),
    zero(dim, zero_f),
    w_grad(dim, boundary_gradw), psi_grad(dim, boundary_gradpsi), 
    grad(&fespace, &fespace_v), rot(dim, rot_f), gradpsi(psi),
    rot_psi_grad(rot, psi_grad), rV_aux(v_aux), rV(rot, zero)
{
    //Create the block offsets
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace.TrueVSize();
    block_true_offsets[2] = fespace.TrueVSize();
    block_true_offsets.PartialSum();
  
    //Initialize the corresponding vectors
    Y.Update(block_true_offsets);
    B.Update(block_true_offsets);
  
    //Define local coefficients
    ConstantCoefficient neg(-1.);
  
    //Dirichlet coefficients
    FunctionCoefficient w_coeff(boundary_w);
    FunctionCoefficient psi_coeff(boundary_psi);
  
    //Rotational coupled coefficients 
    ProductCoefficient neg_w(neg, w_coeff);
    InnerProductCoefficient r_inv_hat_w_grad(r_inv_hat, w_grad);
  
    ScalarVectorProductCoefficient neg_psi_grad(neg, psi_grad);
    InnerProductCoefficient r_inv_hat_psi_grad(r_inv_hat, psi_grad);
  
    //Define essential boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0
  
    ess_bdr_w[0] = 0; ess_bdr_w[1] = 0;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);
  
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace.GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);
  
    //Define grid functions
    w =  new ParGridFunction(&fespace);
    w_aux = new ParGridFunction(&fespace);
    w_aux->ProjectCoefficient(w_coeff);
  
    psi = new ParGridFunction(&fespace);
    psi_aux = new ParGridFunction(&fespace);
    psi_aux->ProjectCoefficient(psi_coeff);
  
    v = new ParGridFunction(&fespace_v);
    v_aux = new ParGridFunction(&fespace_v);
    theta_grad_aux = new ParGridFunction(&fespace);

    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();
  
    g = new ParLinearForm(&fespace);
    g->AddDomainIntegrator(new DomainLFIntegrator(neg_w));
    g->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_psi_grad));
    g->AddDomainIntegrator(new DomainLFGradIntegrator(psi_grad));
    g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_psi_grad), ess_bdr_psi);
    g->Assemble();
    g->ParallelAssemble(B.GetBlock(0));
  
    //Define bilinear forms of the system
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator);
    m->Assemble();
    m->EliminateEssentialBCFromDofs(ess_tdof_list_w);
    m->Finalize();
    M = m->ParallelAssemble();
  
    c = new ParMixedBilinearForm(&fespace, &fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator);
    c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c->Assemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();
  
    this->Update_T(Theta);
}

double r_inv_f(const Vector &x){
    return pow(x(0) + epsilon_r, -1);
}

void r_hat_f(const Vector &x, Vector &y){
    y(0) = 1;
    y(1) = 0;
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 0;
}

void boundary_gradw(const Vector &x, Vector &f){
    f(0) = 0;
    f(1) = 0;
}

//Boundary values for psi
double boundary_psi(const Vector &x){
    return Vel*0.5*pow(x(0), 2);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = Vel*x(0);
    f(1) = 0;
}
