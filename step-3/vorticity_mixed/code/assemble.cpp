#include "header.h"

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

void Artic_sea::assemble_system(){
    //Define local coefficients
    FunctionCoefficient f_coeff(f_rhs);
    ConstantCoefficient g_coeff(0.);
    ConstantCoefficient viscosity(1.);
    FunctionCoefficient eta_coeff(porous_constant);

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
    d->Finalize();
    D = d->ParallelAssemble();

    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(viscosity));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();
    (*M) *= -1.;

    c = new ParMixedBilinearForm(fespace_psi, fespace_w);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(eta_coeff));
    c->Assemble();
    c->Finalize();
    C = c->ParallelAssemble();

    Ct = new TransposeOperator(C);

    //Create the complete bilinear operator:
    //
    //   A = [ M  C^t]
    //       [ C   0 ] 
    A = new BlockOperator(block_true_offsets);
    A->SetBlock(0, 0, D);
    A->SetBlock(0, 1, Ct);
    A->SetBlock(1, 0, C);
    A->SetBlock(1, 1, M);
}

double f_rhs(const Vector &x){                 
    return 2;
}

double porous_constant(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return 1e+6;
    else
        return 0;
}
