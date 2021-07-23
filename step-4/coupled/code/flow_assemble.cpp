#include "header.h"

// unitary r vector
void r_hat(const Vector &x, Vector &y);

double r_inv(const Vector &x);

//Boundary values for w
double boundary_w(const Vector &x);

void boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
double boundary_psi(const Vector &x);

void boundary_gradpsi(const Vector &x, Vector &f);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, ParFiniteElementSpace &fespace_v, int dim, int attributes,  const BlockVector X):
  config(config),
  fespace(fespace),
  block_true_offsets(3),
  f(NULL), g(NULL),
  m(NULL), d(NULL), c(NULL),
  M(NULL), D(NULL), C(NULL), hBlocks(2,2), H(NULL),
  psi(NULL), w(NULL), v(NULL),
  w_aux(NULL), psi_aux(NULL), v_aux(NULL), theta(NULL), phi(NULL),
  r_invCoeff(r_inv), r_hatCoeff(dim, r_hat), r_inv_hat(r_invCoeff, r_hatCoeff),
  neg(-1.), w_coeff(boundary_w), w_grad(dim, boundary_gradw), psi_coeff(boundary_psi),
  psi_grad(dim, boundary_gradpsi), eta_r_inv_hat(neg, r_inv_hat), neg_w(neg, w_coeff),
  r_inv_hat_w_grad(r_inv_hat, w_grad), neg_w_grad(neg, w_grad), r_inv_hat_psi_grad(r_inv_hat, psi_grad),
  eta_r_inv_hat_psi_grad(eta_r_inv_hat, psi_grad), neg_psi_grad(neg, psi_grad),
  eta_psi_grad(neg, psi_grad), neg_eta_psi_grad(neg, eta_psi_grad),
  r(r_f), zero(dim, zero_f), grad(&fespace, &fespace_v), rot(dim, rot_f),
  gradpsi(psi), rot_psi_grad(rot, psi_grad), rV_aux(v_aux), rV(rot, zero),
  ess_bdr_psi(attributes), ess_bdr_w(attributes), superlu(NULL), SLU_A(NULL)
{

  //Create the block offsets
  block_true_offsets[0] = 0;
  block_true_offsets[1] = fespace.TrueVSize();
  block_true_offsets[2] = fespace.TrueVSize();
  block_true_offsets.PartialSum();

  //Initialize the corresponding vectors
  Y.Update(block_true_offsets);
  B.Update(block_true_offsets);

  theta = new ParGridFunction(&fespace);
  theta->SetFromTrueDofs(X.GetBlock(0));

  phi = new ParGridFunction(&fespace);
  phi->SetFromTrueDofs(X.GetBlock(1));

  for (int ii = 0; ii < theta->Size(); ii++){
      (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_fun((*phi)(ii)))));
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

  ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
  ess_bdr_w[2] = 1; ess_bdr_w[3] = 0;
  fespace.GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

  ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
  ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 0;
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

  this->Update_T(X);

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


void Flow_Operator::Update_T(const BlockVector X){
  if(theta) delete theta;
  theta = new ParGridFunction(&fespace);
  theta->SetFromTrueDofs(X.GetBlock(0));

  if(phi) delete phi;
  phi = new ParGridFunction(&fespace);
  phi->SetFromTrueDofs(X.GetBlock(1));

  for (int ii = 0; ii < theta->Size(); ii++){
    (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_fun((*phi)(ii)))));
    (*theta)(ii) = config.epsilon_eta + pow(1-(*theta)(ii), 2)/(pow((*theta)(ii), 3) + config.epsilon_eta);
  }
  T = new ParGridFunction(&fespace);
  T->SetFromTrueDofs(X.GetBlock(0));
  S = new ParGridFunction(&fespace);
  S->SetFromTrueDofs(X.GetBlock(1));

  double a0=8.246121,a1=-4.109344e-2,a2=7.731839e-4,a3=-8.323176e-6,a4=5.516023e-8;
  double b0=-1.827390e-1,b1=3.327135e-3,b2=-5.600674e-5;
  double C=4.886169-2,A,BB;

  Rho = new ParGridFunction(&fespace);
  Rho->SetFromTrueDofs(X.GetBlock(0));
  
  for (int ii = 0; ii < T->Size(); ii++){
    if ((*Rho)(ii) > config.T_f)
      {
	A=a0+a1*(*T)(ii)+a2*pow((*T)(ii),2)+a3*pow((*T)(ii),3)+a4*pow((*T)(ii),4);
	BB=b0+b1*pow((*T)(ii),1)+b2*pow((*T)(ii),2);
	(*Rho)(ii) = A*(*S)(ii)+BB*pow((*S)(ii),1.5)+C*pow((*S)(ii),2);
      }
    else
      (*Rho)(ii) = 0.1;
  }

  //RHS coefficients
  GradientGridFunctionCoefficient grad_theta(T);
  GridFunctionCoefficient grad_rho(Rho);
  InnerProductCoefficient r_hat_grad_theta(r_hatCoeff, grad_theta);
  ProductCoefficient grad_rho_theta(grad_rho, r_hat_grad_theta);
  ProductCoefficient F(-0.31, grad_rho_theta);
  ProductCoefficient rF(r, F);
  ProductCoefficient neg_rF(neg, rF);
  //Define local coefficients
  //Properties coefficients
  eta.SetGridFunction(theta);

  //Rotational coupled coefficients
  eta_r_inv_hat.SetACoef(eta);

  eta_r_inv_hat_psi_grad.SetACoef(eta_r_inv_hat);
  eta_psi_grad.SetACoef(eta);
  neg_eta_psi_grad.SetBCoef(eta_psi_grad);

  if (f) delete f;
  f = new ParLinearForm(&fespace);
  f->AddDomainIntegrator(new DomainLFIntegrator(neg_rF));
  f->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_w_grad));
  f->AddDomainIntegrator(new DomainLFIntegrator(eta_r_inv_hat_psi_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(w_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad), ess_bdr_psi);
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_w_grad), ess_bdr_w);
  f->Assemble();
  f->ParallelAssemble(B.GetBlock(1));

  if (d) delete d;
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
    if (x(0) < entrance)
        return 0.5*vel*pow(x(0), 2)*x(1)/Zmax;
    else
        return 0.5*vel*pow(entrance, 2)*x(1)/Zmax;
}

void boundary_gradpsi(const Vector &x, Vector &f){
    if (x(0) < entrance){
        f(0) = vel*x(0)*x(1)/Zmax;
        f(1) = 0.5*vel*pow(x(0), 2)/Zmax;
    } else {
        f(0) = 0.;
        f(1) = 0.5*vel*pow(entrance, 2)/Zmax;
    }
}
