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
    ProductCoefficient neg_f_coeff(-1., f_coeff);
    ConstantCoefficient viscosity(1.);
    ProductCoefficient neg_viscosity(-1., viscosity);
    FunctionCoefficient eta_coeff(porous_constant);
    ProductCoefficient neg_eta_coeff(-1., eta_coeff);

    //Define grid functions and apply essential boundary conditions(?)
    Array<int> ess_tdof_list_w;
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w = 0; 
    ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 1;
    fespace_w->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

    w =  new ParGridFunction;
    w->MakeRef(fespace_w, x.GetBlock(0), 0);
    w->ProjectBdrCoefficient(boundary_w_coeff, ess_bdr_w);
    w->ParallelProject(x.GetBlock(0));
    x.GetBlock(0).SyncAliasMemory(x);

    Array<int> ess_tdof_list_psi;
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi = 0; 
    ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
    fespace_psi->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    psi = new ParGridFunction;
    psi->MakeRef(fespace_psi, x.GetBlock(1), 0);
    psi->ProjectBdrCoefficient(boundary_psi_coeff, ess_bdr_psi);
    psi->ParallelProject(x.GetBlock(1));
    x.GetBlock(1).SyncAliasMemory(x);

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
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_f_coeff));
    f->Assemble();
    f->SyncAliasMemory(b);
    f->ParallelAssemble(B.GetBlock(1));
    B.GetBlock(1).SyncAliasMemory(B);

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(viscosity));
    m->Assemble();
    m->EliminateEssentialBCFromDofs(ess_tdof_list_w);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace_psi);
    d->AddDomainIntegrator(new DiffusionIntegrator(neg_eta_coeff));
    d->Assemble();
    d->EliminateEssentialBCFromDofs(ess_tdof_list_psi);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace_psi, fespace_w);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(neg_viscosity));
    c->Assemble();
    //c->EliminateTrialDofs(ess_bdr_psi, *psi, *g);
    //c->EliminateTestDofs(ess_bdr_w);
    //c->Finalize();
    //C = c->ParallelAssemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();

    ct = new ParMixedBilinearForm(fespace_w, fespace_psi);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator(neg_viscosity));
    ct->Assemble();
    //ct->EliminateTrialDofs(ess_bdr_w, *w, *f);
    //ct->EliminateTestDofs(ess_bdr_psi);
    //ct->Finalize();
    //Ct = ct->ParallelAssemble();
    OperatorHandle Cth;
    ct->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Cth);
    Ct = Cth.Is<HypreParMatrix>();

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
    return 1.*((out_rad - int_rad)/2 - x(0));
}

double porous_constant(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    double result = 0.;
    if (r_2 < pow(sigma, 2))
        return 1e+6;
    else
        return 0.1;
}
