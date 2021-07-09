#include "header.h"

// unitary r vector
void r_vec(const Vector &x, Vector &y);

//Temperature field
double temperature_f(const Vector &x);

//Right hand side of the equation
 double f_rhs(const Vector &x);

//Boundary values for w
 double scaled_boundary_w(const Vector &x);

 void scaled_boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
 double scaled_boundary_psi(const Vector &x);

 void scaled_boundary_gradpsi(const Vector &x, Vector &f);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, const HypreParVector *X_T):
  fespace(fespace),
  block_true_offsets(3),
  f(NULL), g(NULL),
  m(NULL), d(NULL), c(NULL),
  M(NULL), D(NULL), C(NULL),
  A(NULL), psi(NULL), w(NULL),
  w_aux(NULL), psi_aux(NULL), theta(NULL),
  bg(0.002), f_integrator(NULL), x_T(NULL)
{
  //Create the block offsets
  block_true_offsets[0] = 0;
  block_true_offsets[1] = fespace.TrueVSize();
  block_true_offsets[2] = fespace.TrueVSize();
  block_true_offsets.PartialSum();

  //Initialize the corresponding vectors
  Y.Update(block_true_offsets);
  B.Update(block_true_offsets);

  x_T = new ParGridFunction(&fespace);

  theta = new ParGridFunction(&fespace);
  FunctionCoefficient temperature(temperature_f);
  theta->SetFromTrueDofs(*X_T);
  for (int ii = 0; ii < theta->Size(); ii++){
      (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
      (*theta)(ii) = config.cold_porosity + pow(1-(*theta)(ii), 2)/(pow((*theta)(ii), 3) + config.cold_porosity);
  }
  GridFunctionCoefficient eta(theta);

  //Define local coefficients
  ConstantCoefficient mu(config.viscosity);
  ProductCoefficient neg_mu(-1., mu);
  ProductCoefficient neg_eta(-1., eta);

  FunctionCoefficient w_coeff(scaled_boundary_w);
  VectorFunctionCoefficient w_grad(dim, scaled_boundary_gradw);
  ProductCoefficient neg_mu_w(neg_mu, w_coeff);
  ScalarVectorProductCoefficient mu_w_grad(mu, w_grad);
  ScalarVectorProductCoefficient neg_mu_w_grad(neg_mu, w_grad);

  FunctionCoefficient psi_coeff(scaled_boundary_psi);
  VectorFunctionCoefficient psi_grad(dim, scaled_boundary_gradpsi);
  ScalarVectorProductCoefficient mu_psi_grad(mu, psi_grad);
  ScalarVectorProductCoefficient neg_mu_psi_grad(neg_mu, psi_grad);
  ScalarVectorProductCoefficient eta_psi_grad(eta, psi_grad);
  ScalarVectorProductCoefficient neg_eta_psi_grad(neg_eta, psi_grad);

  //Define essential boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0

  Array<int> ess_tdof_list_w;
  Array<int> ess_bdr_w(attributes);
  ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
  ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
  fespace.GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

  Array<int> ess_tdof_list_psi;
  Array<int> ess_bdr_psi(attributes);
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

  g = new ParLinearForm(&fespace);
  g->AddDomainIntegrator(new DomainLFIntegrator(neg_mu_w));
  g->AddDomainIntegrator(new DomainLFGradIntegrator(mu_psi_grad));
  g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_psi_grad));
  g->Assemble();
  g->ParallelAssemble(B.GetBlock(0));

  //Update the RHS

  //Create the temperature coeffcient (rF)
  x_T->SetFromTrueDofs(*X_T);
  GradientGridFunctionCoefficient delta_T(x_T);
  VectorFunctionCoefficient rcap(2, r_vec);
  InnerProductCoefficient r_deltaT(rcap, delta_T);
  ProductCoefficient bg_deltaT(bg, r_deltaT);
  FunctionCoefficient r(r_f);
  ProductCoefficient rF(bg_deltaT, r);

  f = new ParLinearForm(&fespace);
  f->AddDomainIntegrator(new DomainLFIntegrator(rF));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(mu_w_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_w_grad));
  f->Assemble();
  f->ParallelAssemble(B.GetBlock(1));

  this->Update_T(config, X_T, dim, attributes);

  //Define bilinear forms of the system
  m = new ParBilinearForm(&fespace);
  m->AddDomainIntegrator(new MassIntegrator(mu));
  m->Assemble();
  m->EliminateEssentialBCFromDofs(ess_tdof_list_w, *w, *g);
  m->Finalize();
  M = m->ParallelAssemble();

  d = new ParBilinearForm(&fespace);
  d->AddDomainIntegrator(new DiffusionIntegrator(eta));
  d->Assemble();
  d->EliminateEssentialBCFromDofs(ess_tdof_list_psi, *psi, *f);
  d->Finalize();
  D = d->ParallelAssemble();

  c = new ParMixedBilinearForm(&fespace, &fespace);
  c->AddDomainIntegrator(new MixedGradGradIntegrator(mu));
  c->Assemble();
  OperatorHandle Ch;
  c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
  C = Ch.Is<HypreParMatrix>();
}

void Flow_Operator::Update_T(Config config, const HypreParVector *X_T, int dim, int attributes){

  //Create the temperature coeffcient (rF)
  /*x_T->SetFromTrueDofs(*X_T);
  GradientGridFunctionCoefficient delta_T(x_T);
  VectorFunctionCoefficient rcap(2, r_vec);
  InnerProductCoefficient r_deltaT(rcap, delta_T);
  ProductCoefficient bg_deltaT(bg, r_deltaT);
  FunctionCoefficient r(r_f);
  ProductCoefficient rF(bg_deltaT, r);*/

  //if(f_integrator) delete f_integrator;
  //f_integrator = new DomainLFIntegrator(rF);

  //(f->GetDLFI())[0]=f_integrator;

  //f->Assemble();
  //f->ParallelAssemble(B.GetBlock(1));*/

  if(theta) delete theta;
  theta = new ParGridFunction(&fespace);
  FunctionCoefficient temperature(temperature_f);
  theta->SetFromTrueDofs(*X_T);
  for (int ii = 0; ii < theta->Size(); ii++){
      (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
      (*theta)(ii) = config.cold_porosity + pow(1-(*theta)(ii), 2)/(pow((*theta)(ii), 3) + config.cold_porosity);
  }
  GridFunctionCoefficient eta(theta);

  //Define local coefficients
  ConstantCoefficient mu(config.viscosity);
  ProductCoefficient neg_mu(-1., mu);
  ProductCoefficient neg_eta(-1., eta);

  FunctionCoefficient w_coeff(scaled_boundary_w);
  VectorFunctionCoefficient w_grad(dim, scaled_boundary_gradw);
  ProductCoefficient neg_mu_w(neg_mu, w_coeff);
  ScalarVectorProductCoefficient mu_w_grad(mu, w_grad);
  ScalarVectorProductCoefficient neg_mu_w_grad(neg_mu, w_grad);

  FunctionCoefficient psi_coeff(scaled_boundary_psi);
  VectorFunctionCoefficient psi_grad(dim, scaled_boundary_gradpsi);
  ScalarVectorProductCoefficient mu_psi_grad(mu, psi_grad);
  ScalarVectorProductCoefficient neg_mu_psi_grad(neg_mu, psi_grad);
  ScalarVectorProductCoefficient eta_psi_grad(eta, psi_grad);
  ScalarVectorProductCoefficient neg_eta_psi_grad(neg_eta, psi_grad);

  x_T->SetFromTrueDofs(*X_T);
  GradientGridFunctionCoefficient delta_T(x_T);
  VectorFunctionCoefficient rcap(dim, r_vec);
  InnerProductCoefficient r_deltaT(rcap, delta_T);
  ProductCoefficient bg_deltaT(bg, r_deltaT);
  FunctionCoefficient r(r_f);
  ProductCoefficient rF(bg_deltaT, r);

  if(f) delete f;
  f = new ParLinearForm(&fespace);
  f->AddDomainIntegrator(new DomainLFIntegrator(rF));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
  f->AddDomainIntegrator(new DomainLFGradIntegrator(mu_w_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad));
  f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_w_grad));
  f->Assemble();
  f->ParallelAssemble(B.GetBlock(1));

  Array<int> ess_tdof_list_psi;
  Array<int> ess_bdr_psi(attributes);
  ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
  ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
  fespace.GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

  if(d) delete d;
  d = new ParBilinearForm(&fespace);
  d->AddDomainIntegrator(new DiffusionIntegrator(eta));
  d->Assemble();
  d->EliminateEssentialBCFromDofs(ess_tdof_list_psi, *psi, *f);
  d->Finalize();
  D = d->ParallelAssemble();
}


void r_vec(const Vector &x, Vector &y){
  y(0)=1;
  y(1)=0;
}

//Temperature field
double temperature_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = height/2;
    double sigma = (Rmax - Rmin)/10;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;

}

//Right hand side of the equation
double f_rhs(const Vector &x){
    return 0.;
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 0.;
}

void boundary_gradw(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

//Boundary values for psi
double boundary_psi(const Vector &x){
    return x(0);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = 1.;
    f(1) = 0.;
}

//Scaling for the boundary conditions
double left_border(const Vector &x){
    double width = Rmax - Rmin;
    return 0.5*(1 - tanh((5*border/width)*(Rmin + width/border - x(0))));
}

double right_border(const Vector &x){
    double width = Rmax - Rmin;
    return 0.5*(1 - tanh((5*border/width)*(x(0) - Rmax + width/border)));
}

double lower_border(const Vector &x){
    return 0.5*(1 - tanh((5*border/height)*(height/border - x(1))));
}

double upper_border(const Vector &x){
    return 0.5*(1 - tanh((5*border/height)*(x(1) - height + height/border)));
}

double left_grad(const Vector &x){
    double width = Rmax - Rmin;
    return (2.5*border/width)*(1 - pow(tanh((5*border/width)*(Rmin + width/border - x(0))), 2));
}

double right_grad(const Vector &x){
    double width = Rmax - Rmin;
    return -(2.5*border/width)*(1 - pow(tanh((5*border/width)*(x(0) - Rmax + width/border)), 2));
}

double lower_grad(const Vector &x){
    return (2.5*border/height)*(1 - pow(tanh((5*border/height)*(height/border - x(1))), 2));
}

double upper_grad(const Vector &x){
    return -(2.5*border/height)*(1 - pow(tanh((5*border/height)*(x(1) - height + height/border)), 2));
}

double scaled_boundary_w(const Vector &x){
    return (1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x))*boundary_w(x);
}

void scaled_boundary_gradw(const Vector &x, Vector &f){
    boundary_gradw(x, f);

    f(0) *= 1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x);
    f(1) *= 1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x);

    f(0) -= boundary_w(x)*(left_grad(x)*right_border(x) + left_border(x)*right_grad(x))*lower_border(x)*upper_border(x);
    f(1) -= boundary_w(x)*left_border(x)*right_border(x)*(lower_grad(x)*upper_border(x) + lower_border(x)*upper_grad(x));
}

double scaled_boundary_psi(const Vector &x){
    return (1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x))*boundary_psi(x);
}

void scaled_boundary_gradpsi(const Vector &x, Vector &f){
    boundary_gradpsi(x, f);

    f(0) *= 1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x);
    f(1) *= 1 - left_border(x)*right_border(x)*lower_border(x)*upper_border(x);

    f(0) -= boundary_psi(x)*(left_grad(x)*right_border(x) + left_border(x)*right_grad(x))*lower_border(x)*upper_border(x);
    f(1) -= boundary_psi(x)*left_border(x)*right_border(x)*(lower_grad(x)*upper_border(x) + lower_border(x)*upper_grad(x));
}
