#include "header.h"

//Rotational functions
void r_inv_hat_f(const Vector &x, Vector &f);

//Temperature field
double temperature_f(const Vector &x);

//Right hand side of the equation
double f_rhs(const Vector &x);

//Boundary values for w
double boundary_w(const Vector &x);

void boundary_gradw(const Vector &x, Vector &f);

//Boundary values for psi
double boundary_psi(const Vector &x);

void boundary_gradpsi(const Vector &x, Vector &f);

void Artic_sea::assemble_system(){
    //Calculate the porus coefficient
    theta = new ParGridFunction(fespace);
    FunctionCoefficient temperature(temperature_f);
    theta->ProjectCoefficient(temperature);
    for (int ii = 0; ii < theta->Size(); ii++){
        (*theta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - config.T_f)));
        (*theta)(ii) = config.epsilon_eta + (1 - pow((*theta)(ii), 2))/(pow((*theta)(ii), 3) + config.epsilon_eta);
    }

    //Rotational coefficients
    VectorFunctionCoefficient r_inv_hat(dim, r_inv_hat_f);

    //Properties coefficients
    ConstantCoefficient inv_mu(pow(config.viscosity, -1));
    GridFunctionCoefficient eta(theta);
    ConstantCoefficient neg(-1.);

    //RHS coefficients
    FunctionCoefficient f_coeff(f_rhs);
    ProductCoefficient inv_mu_f(inv_mu, f_coeff);

    //Dirichlet coefficients
    FunctionCoefficient w_coeff(boundary_w);
    VectorFunctionCoefficient w_grad(dim, boundary_gradw);

    FunctionCoefficient psi_coeff(boundary_psi);
    VectorFunctionCoefficient psi_grad(dim, boundary_gradpsi);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient eta_r_inv_hat(eta, r_inv_hat);

    ProductCoefficient neg_w(neg, w_coeff);
    InnerProductCoefficient r_inv_hat_w_grad(r_inv_hat, w_grad);
    ScalarVectorProductCoefficient neg_w_grad(neg, w_grad);

    InnerProductCoefficient r_inv_hat_psi_grad(r_inv_hat, psi_grad);
    InnerProductCoefficient eta_r_inv_hat_psi_grad(eta_r_inv_hat, psi_grad);
    ScalarVectorProductCoefficient neg_psi_grad(neg, psi_grad);
    ScalarVectorProductCoefficient eta_psi_grad(eta, psi_grad);
    ScalarVectorProductCoefficient neg_eta_psi_grad(neg, eta_psi_grad);

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
    ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

    Array<int> ess_tdof_list_psi;
    ess_bdr_psi.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    //Define grid functions
    w =  new ParGridFunction(fespace);
    w_aux = new ParGridFunction(fespace);
    w_aux->ProjectCoefficient(w_coeff);

    psi = new ParGridFunction(fespace);
    psi_aux = new ParGridFunction(fespace);
    psi_aux->ProjectCoefficient(psi_coeff);

    //Define the RHS
    g = new ParLinearForm(fespace);
    g->AddDomainIntegrator(new DomainLFIntegrator(neg_w));
    g->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_psi_grad));
    g->AddDomainIntegrator(new DomainLFGradIntegrator(psi_grad));
    g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_psi_grad), ess_bdr_psi);
    g->Assemble();
    g->ParallelAssemble(B.GetBlock(0));

    f = new ParLinearForm(fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(inv_mu_f));
    f->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_w_grad));
    f->AddDomainIntegrator(new DomainLFIntegrator(eta_r_inv_hat_psi_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(w_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_w_grad), ess_bdr_w);
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad), ess_bdr_psi);
    f->Assemble();
    f->ParallelAssemble(B.GetBlock(1));

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator);
    m->Assemble();
    m->EliminateEssentialBCFromDofs(ess_tdof_list_w);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta));
    d->AddDomainIntegrator(new ConvectionIntegrator(eta_r_inv_hat));
    d->Assemble();
    d->EliminateEssentialBCFromDofs(ess_tdof_list_psi);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace, fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator);
    c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c->Assemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();
}

double r_f(const Vector &x){
    return x(0);
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0) + epsilon_r, -1);
    f(1) = 0.;
}

//Temperature field
double temperature_f(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/10;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -0.01;
    else
        return 0;
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
    return 0.5*pow(x(0), 2);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = x(0);
    f(1) = 0.;
}
