#include "header.h"

//Boundary values for w
double boundary_w(const Vector &x);
void boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
double boundary_psi(const Vector &x);
void boundary_gradpsi(const Vector &x, Vector &f);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

void Artic_sea::assemble_system(){
    double viscosity = 1.;

    //Define local coefficients
    ConstantCoefficient zero(0.);
    ConstantCoefficient mu(viscosity);
    ProductCoefficient neg_mu(-1., mu);
    FunctionCoefficient eta(porous_constant);
    ProductCoefficient neg_eta(-1., eta);

    FunctionCoefficient f_coeff(f_rhs);
    ProductCoefficient neg_f_coeff(-1., f_coeff);
    FunctionCoefficient w_aux(boundary_w);
    VectorFunctionCoefficient gradw_aux(dim, boundary_gradw);
    FunctionCoefficient psi_aux(boundary_psi);
    VectorFunctionCoefficient gradpsi_aux(dim, boundary_gradpsi);

    ProductCoefficient neg_mu_w_aux(neg_mu, w_aux);
    ScalarVectorProductCoefficient mu_gradpsi_aux(mu, gradpsi_aux);
    ScalarVectorProductCoefficient neg_mu_gradpsi_aux(neg_mu, gradpsi_aux);
    ScalarVectorProductCoefficient mu_gradw_aux(mu, gradw_aux);
    ScalarVectorProductCoefficient neg_mu_gradw_aux(neg_mu, gradw_aux);
    ScalarVectorProductCoefficient eta_gradpsi_aux(eta, gradpsi_aux);
    ScalarVectorProductCoefficient neg_eta_gradpsi_aux(neg_eta, gradpsi_aux);
    
    ConstantCoefficient boundary_w_coeff(0.);
    FunctionCoefficient boundary_psi_coeff(boundary_psi);
    ConstantCoefficient g_coeff(0.);

    //Define grid functions
    w =  new ParGridFunction(fespace_w);
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w = 0; 
    ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 1;
    w->ProjectBdrCoefficient(zero, ess_bdr_w);
    w->ParallelProject(x.GetBlock(0));

    psi = new ParGridFunction(fespace_psi);
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi = 0; 
    ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
    psi->ProjectBdrCoefficient(zero, ess_bdr_psi);
    psi->ParallelProject(x.GetBlock(1));

    //Define essential boundary conditions
    w_boundary = new ParGridFunction(fespace_w);
    w_boundary->ProjectCoefficient(w_aux);

    psi_boundary = new ParGridFunction(fespace_psi);
    psi_boundary->ProjectCoefficient(psi_aux);

    //Define the RHS
    g = new ParLinearForm;
    g->Update(fespace_w, b.GetBlock(0), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(neg_mu_w_aux));
    g->AddDomainIntegrator(new DomainLFGradIntegrator(mu_gradpsi_aux));
    g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_gradpsi_aux));
    g->Assemble();
    g->SyncAliasMemory(b);
    g->ParallelAssemble(B.GetBlock(0));
    B.GetBlock(0).SyncAliasMemory(B);

    f = new ParLinearForm;
    f->Update(fespace_psi, b.GetBlock(1), 0);
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_f_coeff));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(mu_gradw_aux));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_gradpsi_aux));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_gradw_aux));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_gradpsi_aux));
    f->Assemble();
    f->SyncAliasMemory(b);
    f->ParallelAssemble(B.GetBlock(1));
    B.GetBlock(1).SyncAliasMemory(B);

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(mu));
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace_psi);
    d->AddDomainIntegrator(new DiffusionIntegrator(neg_eta));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace_psi, fespace_w);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(neg_mu));
    c->Assemble();
    c->EliminateTrialDofs(ess_bdr_psi, *psi, *f);
    c->EliminateTestDofs(ess_bdr_w);
    c->Finalize();
    C = c->ParallelAssemble();

    ct = new ParMixedBilinearForm(fespace_w, fespace_psi);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator(neg_mu));
    ct->Assemble();
    ct->EliminateTrialDofs(ess_bdr_w, *w, *g);
    ct->EliminateTestDofs(ess_bdr_psi);
    ct->Finalize();
    Ct = ct->ParallelAssemble();

    //Create the complete bilinear operator:
    //
    //   A = [ M    C ] 
    //       [ C^t  D ] 
    A = new BlockOperator(block_true_offsets);
    A->SetBlock(0, 0, M);
    A->SetBlock(0, 1, C);
    A->SetBlock(1, 0, Ct);
    A->SetBlock(1, 1, D);
}

double boundary_w(const Vector &x){
    return 0.;
}

void boundary_gradw(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

double boundary_psi(const Vector &x){
    return x(0);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = 1.;
    f(1) = 0.;
}

double f_rhs(const Vector &x){                 
    return 0;
}

double porous_constant(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    double result = 0.;
    if (r_2 < pow(sigma, 2))
        result = 1e+12;
    else
        result = 0.1;
    return -result;
}
