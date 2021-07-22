#include "header.h"

// unitary r vector
void r_hat(const Vector &x, Vector &y);

//void r_inv_hat_f(const Vector &x, Vector &f);

double r_inv(const Vector &x);

//Right hand side of the equation
 double f_rhs(const Vector &x);

//Boundary values for w
 double boundary_w(const Vector &x);

 void boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
 double boundary_psi(const Vector &x);

 void boundary_gradpsi(const Vector &x, Vector &f);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes, Array<int> block_true_offsets, const HypreParVector *X_T):
  fespace(fespace),
  f(NULL), g(NULL),
  m(NULL), d(NULL), c(NULL),
  M(NULL), D(NULL), C(NULL), hBlocks(2,2), H(NULL),
  psi(NULL), w(NULL), v(NULL),
  w_aux(NULL), psi_aux(NULL), theta(NULL),
  bg(1.), r_invCoeff(r_inv), r_hatCoeff(dim, r_hat), r_inv_hat(r_invCoeff, r_hatCoeff),
  neg(-1.), w_coeff(boundary_w), w_grad(dim, boundary_gradw), psi_coeff(boundary_psi),
  psi_grad(dim, boundary_gradpsi), eta_r_inv_hat(neg, r_inv_hat), neg_w(neg, w_coeff),
  r_inv_hat_w_grad(r_inv_hat, w_grad), neg_w_grad(neg, w_grad), r_inv_hat_psi_grad(r_inv_hat, psi_grad),
  eta_r_inv_hat_psi_grad(eta_r_inv_hat, psi_grad), neg_psi_grad(neg, psi_grad),
  eta_psi_grad(neg, psi_grad), neg_eta_psi_grad(neg, eta_psi_grad), inv_mu(pow(config.viscosity, -1)),
  x_T(NULL), delta_T(x_T), rcap(dim, r_hat), r_deltaT(rcap, delta_T), bg_deltaT(bg, r_deltaT),
  r(r_f), rF(bg_deltaT, r), inv_mu_TCoeff(inv_mu, rF),
  ess_bdr_psi(attributes), ess_bdr_w(attributes), superlu(NULL), SLU_A(NULL)
{
  //Initialize the corresponding vectors
  Y.Update(block_true_offsets);
  B.Update(block_true_offsets);

  x_T = new ParGridFunction(&fespace);

  theta = new ParGridFunction(&fespace);
  theta->SetFromTrueDofs(*X_T);
  for (int ii = 0; ii < theta->Size(); ii++){
      (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
      (*theta)(ii) = config.epsilon_eta + pow(1-(*theta)(ii), 2)/(pow((*theta)(ii), 3) + config.epsilon_eta);
  }

  //Define local coefficients
  //Properties coefficients
  eta.SetGridFunction(theta);

  //Rotational coupled coefficients
  eta_r_inv_hat.SetACoef(eta);

  eta_r_inv_hat_psi_grad.SetACoef(eta_r_inv_hat);
  eta_psi_grad.SetACoef(eta);
  neg_eta_psi_grad.SetBCoef(eta_psi_grad);

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

  ess_bdr_psi[0] = 0; ess_bdr_psi[1] = 0;
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

  g = new ParLinearForm(&fespace);
  g->AddDomainIntegrator(new DomainLFIntegrator(neg_w));
  g->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_psi_grad));
  g->AddDomainIntegrator(new DomainLFGradIntegrator(psi_grad));
  g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_psi_grad), ess_bdr_psi);
  g->Assemble();
  g->ParallelAssemble(B.GetBlock(0));

  this->Update_T(config, X_T, dim, attributes);

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
}


void Flow_Operator::Update_T(Config config, const HypreParVector *X_T, int dim, int attributes){
  if(theta) delete theta;
  theta = new ParGridFunction(&fespace);
  theta->SetFromTrueDofs(*X_T);
  for (int ii = 0; ii < theta->Size(); ii++){
    (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
    (*theta)(ii) = config.epsilon_eta + pow(1-(*theta)(ii), 2)/(pow((*theta)(ii), 3) + config.epsilon_eta);
  }

  //Define local coefficients
  //Properties coefficients
  eta.SetGridFunction(theta);

  //Rotational coupled coefficients
  eta_r_inv_hat.SetACoef(eta);

  eta_r_inv_hat_psi_grad.SetACoef(eta_r_inv_hat);
  eta_psi_grad.SetACoef(eta);
  neg_eta_psi_grad.SetBCoef(eta_psi_grad);

  x_T->SetFromTrueDofs(*X_T);
  delta_T.SetGridFunction(x_T);
  r_deltaT.SetBCoef(delta_T);
  bg_deltaT.SetBCoef(r_deltaT);
  rF.SetACoef(bg_deltaT);
  inv_mu_TCoeff.SetBCoef(rF);


  if(f) delete f;
  f = new ParLinearForm(&fespace);
  f->AddDomainIntegrator(new DomainLFIntegrator(inv_mu_TCoeff));
  f->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_w_grad));
  f->AddDomainIntegrator(new DomainLFIntegrator(eta_r_inv_hat_psi_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(w_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad), ess_bdr_psi);
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_w_grad), ess_bdr_w);
  f->Assemble();
  f->ParallelAssemble(B.GetBlock(1));

  if(d) delete d;
  d = new ParBilinearForm(&fespace);
  d->AddDomainIntegrator(new DiffusionIntegrator(eta));
  d->AddDomainIntegrator(new ConvectionIntegrator(eta_r_inv_hat));
  d->Assemble();
  d->EliminateEssentialBCFromDofs(ess_tdof_list_psi);
  d->Finalize();
  D = d->ParallelAssemble();
}

double r_inv(const Vector &x){
  return pow(x(0) + epsilon_r, -1);
}

void r_hat(const Vector &x, Vector &y){
  y(0)=1;
  y(1)=0;
}

//Right hand side of the equation
double f_rhs(const Vector &x){
    return 0.;
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
    return 0.5*x(0)*x(0);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = x(0);
    f(1) = 0;
}
