#include "header.h"

//Right hand side of the equation
double f_rhs(const Vector &x);

//Constant to the brinkman term
double porous_constant(const Vector &x);

void Artic_sea::assemble_system(){
    //Define local coefficients
    ConstantCoefficient g_coeff(0.);
    FunctionCoefficient f_coeff(f_rhs);

    ConstantCoefficient viscosity(1.);
    FunctionCoefficient eta_coeff(porous_constant);

    ProductCoefficient neg_f_coeff(-1., f_coeff);

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
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w[0] = 0; ess_bdr_w[1] = 0;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
    fespace_w->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

    Array<int> ess_tdof_list_psi;
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 0; ess_bdr_psi[1] = 0;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace_psi->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    //Define grid functions
    w =  new ParGridFunction;
    w->MakeRef(fespace_w, x.GetBlock(0), 0);

    psi = new ParGridFunction;
    psi->MakeRef(fespace_psi, x.GetBlock(1), 0);

    //Define the RHS
    g = new ParLinearForm;
    g->Update(fespace_w, b.GetBlock(0), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(g_coeff));
    g->Assemble();
    g->ParallelAssemble(B.GetBlock(0));

    f = new ParLinearForm;
    f->Update(fespace_psi, b.GetBlock(1), 0);
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_f_coeff));
    f->Assemble();
    f->ParallelAssemble(B.GetBlock(1));

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace_w);
    m->AddDomainIntegrator(new MassIntegrator(viscosity));
    m->Assemble();
    m->EliminateEssentialBCFromDofs(ess_tdof_list_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace_psi);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta_coeff));
    d->Assemble();
    d->EliminateEssentialBCFromDofs(ess_tdof_list_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace_psi, fespace_w);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(viscosity));
    c->Assemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();
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
