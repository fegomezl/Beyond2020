#include "header.h"

//Boundary values for psi
double boundary_psi(const Vector &x);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

//r as a vector function
void r_vec(const Vector &x, Vector &y);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace, int attributes, const HypreParVector *X_T):
  fespace(fespace),
  block_offsets(3), block_true_offsets(3),
  f(NULL), g(NULL),
  m(NULL), d(NULL), c(NULL), ct(NULL),
  M(NULL), D(NULL), C(NULL), Ct(NULL),
  A(NULL), psi(NULL), w(NULL), bg(0.002)
{
  //Create the block offsets
  block_offsets[0] = 0;
  block_offsets[1] = fespace.GetVSize();
  block_offsets[2] = fespace.GetVSize();
  block_offsets.PartialSum();

  block_true_offsets[0] = 0;
  block_true_offsets[1] = fespace.TrueVSize();
  block_true_offsets[2] = fespace.TrueVSize();
  block_true_offsets.PartialSum();

  //Initialize the corresponding vectors
  /*const char *device_config = "cpu";
  Device device(device_config);
  MemoryType mt = device.GetMemoryType();*/
  y.Update(block_offsets); Y.Update(block_true_offsets);
  b.Update(block_offsets); B.Update(block_true_offsets);

  FunctionCoefficient boundary_psi_coeff(boundary_psi);
  ConstantCoefficient g_coeff(0.);
  ConstantCoefficient viscosity(1.);
  FunctionCoefficient eta_coeff(porous_constant);

  //Define grid functions and apply essential boundary conditions(?)
  psi = new ParGridFunction(&fespace);
  Array<int> ess_bdr_psi(attributes);
  ess_bdr_psi = 0;
  ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
  psi->ProjectBdrCoefficient(boundary_psi_coeff, ess_bdr_psi);
  psi->ParallelProject(y.GetBlock(0));

  w =  new ParGridFunction(&fespace);
  Array<int> ess_bdr_w(attributes);
  ess_bdr_w = 0;
  ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 1;
  w->ProjectBdrCoefficient(g_coeff, ess_bdr_w);
  w->ParallelProject(y.GetBlock(1));

  this->Update_T(config,X_T);

  g = new ParLinearForm;
  g->Update(&fespace, b.GetBlock(1), 0);
  g->AddDomainIntegrator(new DomainLFIntegrator(g_coeff));
  g->Assemble();
  g->SyncAliasMemory(b);
  g->ParallelAssemble(B.GetBlock(1));
  B.GetBlock(1).SyncAliasMemory(B);

  //Define bilinear forms of the system
  d = new ParBilinearForm(&fespace);
  d->AddDomainIntegrator(new DiffusionIntegrator(eta_coeff));
  d->Assemble();
  d->EliminateEssentialBC(ess_bdr_psi, *psi, *f);
  d->Finalize();
  D = d->ParallelAssemble();

  m = new ParBilinearForm(&fespace);
  m->AddDomainIntegrator(new MassIntegrator(viscosity));
  m->Assemble();
  m->EliminateEssentialBC(ess_bdr_w, *w, *g);
  m->Finalize();
  M = m->ParallelAssemble();
  (*M) *= -1.;

  c = new ParMixedBilinearForm(&fespace, &fespace);
  c->AddDomainIntegrator(new MixedGradGradIntegrator(viscosity));
  c->Assemble();
  c->EliminateTrialDofs(ess_bdr_psi, *psi, *f);
  c->EliminateTestDofs(ess_bdr_w);
  c->Finalize();
  C = c->ParallelAssemble();

  ct = new ParMixedBilinearForm(&fespace, &fespace);
  ct->AddDomainIntegrator(new MixedGradGradIntegrator(viscosity));
  ct->Assemble();
  ct->EliminateTrialDofs(ess_bdr_w, *w, *g);
  ct->EliminateTestDofs(ess_bdr_psi);
  ct->Finalize();
  Ct = ct->ParallelAssemble();
  //Ct = new TransposeOperator(C);

  //Create the complete bilinear operator:
  //
  //   A = [ D  C^t]
  //       [ C  M ]

  A = new BlockOperator(block_true_offsets);
  A->SetBlock(0, 0, D);
  A->SetBlock(0, 1, Ct);
  A->SetBlock(1, 0, C);
  A->SetBlock(1, 1, M);
}

double boundary_psi(const Vector &x){
    return x(0);
}

double f_rhs(const Vector &x){
    return (Rmax - Rmin)/2 - x(0);
}

double porous_constant(const Vector &x){
    double height = Zmax - Zmin;
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = height/2;
    double sigma = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return 1e+6;
    else
        return 0.1;
}

void r_vec(const Vector &x, Vector &y){
  y(0)=1;
  y(1)=0;
}

void Flow_Operator::Update_T(Config config, const HypreParVector *X_T){

  //Create the temperature coeffcient (rF)
  /*GradientGridFunctionCoefficient delta_T(x_T);
  VectorFunctionCoefficient rcap(2, r_vec);
  InnerProductCoefficient r_deltaT(rcap, delta_T);
  ProductCoefficient bg_deltaT(bg, r_deltaT);
  FunctionCoefficient r(r_f);
  ProductCoefficient rF(bg_deltaT, r);*/

  FunctionCoefficient f_coeff(f_rhs);

  //Update the RHS
  if(f) delete f;
  f = new ParLinearForm;
  f->Update(&fespace, b.GetBlock(0), 0);
  f->AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
  f->Assemble();
  f->SyncAliasMemory(b);
  f->ParallelAssemble(B.GetBlock(0));
  B.GetBlock(0).SyncAliasMemory(B);
}
