#include "header.h"

double r_f(const Vector &x);

void r_inv_hat_f(const Vector &x, Vector &f);

//Temperature field
double temperature_f(const Vector &x);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Boundary values for w
double boundary_w(const Vector &x);

//Boundary values for psi
double boundary_psi(const Vector &x);

void Artic_sea::assemble_system(){
    //Calculate the porus coefficient
    theta = new ParGridFunction(fespace);
    FunctionCoefficient temperature(temperature_f);
    theta->ProjectCoefficient(temperature);
    for (int ii = 0; ii < theta->Size(); ii++){
        (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
        (*theta)(ii) = config.cold_porosity + (1 - pow((*theta)(ii), 2))/(pow((*theta)(ii), 3) + config.cold_porosity);
    }

    //Rotational coefficients
    FunctionCoefficient r_coeff(r_f);
    VectorFunctionCoefficient r_inv_hat(dim, r_inv_hat_f);

    //Properties coefficients
    ConstantCoefficient mu(config.viscosity);
    GridFunctionCoefficient eta(theta);
    ProductCoefficient neg_mu(-1., mu);
    ProductCoefficient neg_eta(-1., eta);

    //RHS coefficients
    FunctionCoefficient f_coeff(f_rhs);
    ProductCoefficient neg_f_coeff(-1., f_coeff);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient mu_r_inv_hat(mu, r_inv_hat);
    ScalarVectorProductCoefficient eta_r_inv_hat(eta, r_inv_hat);
    ProductCoefficient neg_r_f_coeff(r_coeff, neg_f_coeff);

    //Dirichlet coefficients
    FunctionCoefficient w_coeff(boundary_w);

    FunctionCoefficient psi_coeff(boundary_psi);

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
    fespace->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

    Array<int> ess_tdof_list_psi;
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 0; ess_bdr_psi[1] = 0;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    //Define grid functions
    w =  new ParGridFunction(fespace);
    w->Randomize();
    w->ProjectBdrCoefficient(w_coeff, ess_bdr_w);
    w->GetTrueDofs(X.GetBlock(0));

    psi = new ParGridFunction(fespace);
    psi->Randomize();
    psi->ProjectBdrCoefficient(psi_coeff, ess_bdr_psi);
    psi->GetTrueDofs(X.GetBlock(1));

    //Define the RHS
    g = new ParLinearForm(fespace);
    g->Assemble();

    f = new ParLinearForm(fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_r_f_coeff));
    f->Assemble();

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(mu));
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta));
    //d->AddDomainIntegrator(new ConvectionIntegrator(eta_r_inv_hat));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace, fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(mu));
    //c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(mu_r_inv_hat));
    c->Assemble();
    c->EliminateTrialDofs(ess_bdr_psi, *psi, *g);
    c->EliminateTestDofs(ess_bdr_w);
    c->Finalize();
    C = c->ParallelAssemble();

    ct = new ParMixedBilinearForm(fespace, fespace);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator(mu));
    //ct->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(mu_r_inv_hat));
    ct->Assemble();
    ct->EliminateTrialDofs(ess_bdr_w, *w, *f);
    ct->EliminateTestDofs(ess_bdr_psi);
    ct->Finalize();
    Ct = ct->ParallelAssemble();

    g->ParallelAssemble(B.GetBlock(0));
    f->ParallelAssemble(B.GetBlock(1));
}

double r_f(const Vector &x){
    return x(0);
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0) + InvR, -1);
    f(1) = 0.;
}

//Temperature field
double temperature_f(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/10;

    //double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;
}

//Right hand side of the equation
double f_rhs(const Vector &x){                 
    return -1.;
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 0.;
}

double vel = 1e-2;

//Boundary values for psi
double boundary_psi(const Vector &x){
    return -vel*0.5*pow(x(0), 2);
    //return (vel*0.5/out_rad)*x(0)*(x(0) - 2*out_rad);
}
