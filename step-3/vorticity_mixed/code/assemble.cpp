#include "header.h"

//Boundary values for psi
double boundary_psi(const Vector &x);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

void Artic_sea::assemble_system(){
    //Define local coefficients
    FunctionCoefficient boundary_psi_coeff(boundary_psi);
    FunctionCoefficient f_coeff(f_rhs);
    ConstantCoefficient g_coeff(0.0);
    ConstantCoefficient viscosity(1.0);
    FunctionCoefficient eta_coeff(porous_constant);

    //Define grid functions and apply essential boundary conditions(?)
    psi = new ParGridFunction(fespace_psi);
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi = 0; 
    ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
    psi->ProjectBdrCoefficient(boundary_psi_coeff, ess_bdr_psi);
    psi->ParallelProject(x.GetBlock(0));

    w =  new ParGridFunction(fespace_w);
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w = 0; 
    ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 1;
    w->ProjectBdrCoefficient(g_coeff, ess_bdr_w);
    w->ParallelProject(x.GetBlock(1));

    //Define the RHS
    f = new ParLinearForm;
    f->Update(fespace_psi, b.GetBlock(0), 0);
    f->AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    f->Assemble();
    f->SyncAliasMemory(b);
    f->ParallelAssemble(B.GetBlock(0));
    B.GetBlock(0).SyncAliasMemory(B);

    g = new ParLinearForm;
    g->Update(fespace_w, b.GetBlock(1), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(g_coeff));
    g->Assemble();
    g->SyncAliasMemory(b);
    g->ParallelAssemble(B.GetBlock(1));
    B.GetBlock(1).SyncAliasMemory(B);

    //Define bilinear forms of the system
    d = new ParBilinearForm(fespace_psi);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta_coeff));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(viscosity));
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();
    (*M) *= -1.;

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
    return (out_rad - int_rad)/2 - x(0);
}

double porous_constant(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return 1e+6;
    else
        return 0.1;
}
