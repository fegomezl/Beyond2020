#include "header.h"

//Boundary values for psi
double boundary_psi(const Vector &x);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

void Artic_sea::assemble_system(){
    //Define local coefficients
    ConstantCoefficient boundary_w_coeff(0.);
    FunctionCoefficient boundary_psi_coeff(boundary_psi);
    ConstantCoefficient g_coeff(0.);
    FunctionCoefficient f_coeff(f_rhs);
    ConstantCoefficient viscosity(-1.);
    FunctionCoefficient eta_coeff(porous_constant);

    //Define grid functions and apply essential boundary conditions(?)
    w =  new ParGridFunction(fespace_w);
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w = 0; 
    ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 20;
    w->ProjectBdrCoefficient(g_coeff, ess_bdr_w);
    w->ParallelProject(x.GetBlock(0));

    psi = new ParGridFunction(fespace_psi);
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi = 0; 
    ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
    psi->ProjectBdrCoefficient(boundary_psi_coeff, ess_bdr_psi);
    psi->ParallelProject(x.GetBlock(1));

    //Define the RHS
    g = new ParLinearForm;
    g->Update(fespace_w, b.GetBlock(0), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(g_coeff));
    g->Assemble();
    g->SyncAliasMemory(b);
    g->ParallelAssemble(B.GetBlock(0));
    B.GetBlock(0).SyncAliasMemory(B);

    f = new ParLinearForm;
    f->Update(fespace_psi, b.GetBlock(1), 0);
    f->AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    f->Assemble();
    f->SyncAliasMemory(b);
    f->ParallelAssemble(B.GetBlock(1));
    B.GetBlock(1).SyncAliasMemory(B);

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(viscosity));
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();
    (*M) *= -1.;

    d = new ParBilinearForm(fespace_psi);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta_coeff));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace_psi, fespace_w);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(viscosity));
    c->Assemble();
    c->EliminateTrialDofs(ess_bdr_psi, *psi, *f);
    c->EliminateTestDofs(ess_bdr_w);
    c->Finalize();
    C = c->ParallelAssemble();

    ct = new ParMixedBilinearForm(fespace_w, fespace_psi);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator(viscosity));
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

double boundary_psi(const Vector &x){
    return x(0);
}

double f_rhs(const Vector &x){                 
    double result = (out_rad - int_rad)/2 - x(0);
    return -result;
}

double porous_constant(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    double result = 0.;
    if (r_2 < pow(sigma, 2))
        result = 1e+6;
    else
        result = 0.1;
    return -result;
}
