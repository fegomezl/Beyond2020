#include "header.h"

//Rotational functions
double r_f(const Vector &x);

void r_inv_hat_f(const Vector &x, Vector &f);

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

    //Dirichlet coefficients
    FunctionCoefficient w_coeff(scaled_boundary_w);
    VectorFunctionCoefficient w_grad(dim, scaled_boundary_gradw);

    FunctionCoefficient psi_coeff(scaled_boundary_psi);
    VectorFunctionCoefficient psi_grad(dim, scaled_boundary_gradpsi);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient mu_r_inv_hat(mu, r_inv_hat);
    ScalarVectorProductCoefficient eta_r_inv_hat(eta, r_inv_hat);

    ProductCoefficient neg_mu_w(neg_mu, w_coeff);
    InnerProductCoefficient mu_r_inv_hat_w_grad(mu_r_inv_hat, w_grad);
    ScalarVectorProductCoefficient mu_w_grad(mu, w_grad);
    ScalarVectorProductCoefficient neg_mu_w_grad(neg_mu, w_grad);

    InnerProductCoefficient mu_r_inv_hat_psi_grad(mu_r_inv_hat, psi_grad);
    InnerProductCoefficient eta_r_inv_hat_psi_grad(eta_r_inv_hat, psi_grad);
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
    g->AddDomainIntegrator(new DomainLFIntegrator(neg_mu_w));
    g->AddDomainIntegrator(new DomainLFIntegrator(mu_r_inv_hat_psi_grad));
    g->AddDomainIntegrator(new DomainLFGradIntegrator(mu_psi_grad));
    g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_psi_grad), ess_bdr_psi);
    g->Assemble();
    g->ParallelAssemble(B.GetBlock(0));

    f = new ParLinearForm(fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    f->AddDomainIntegrator(new DomainLFIntegrator(mu_r_inv_hat_w_grad));
    f->AddDomainIntegrator(new DomainLFIntegrator(eta_r_inv_hat_psi_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(mu_w_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_w_grad), ess_bdr_w);
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad), ess_bdr_psi);
    f->Assemble();
    f->ParallelAssemble(B.GetBlock(1));

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(mu));
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
    c->AddDomainIntegrator(new MixedGradGradIntegrator(mu));
    c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(mu_r_inv_hat));
    c->Assemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();
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
    /*if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;*/
    return 10.;
}

double scale = 1e-5;

//Right hand side of the equation
double f_rhs(const Vector &x){                 
    return 8*scale*pow(x(0), 2)*x(1);
}

//Boundary values for w
double boundary_w(const Vector &x){
    return -8*scale*pow(x(0), 2)*x(1);
}

void boundary_gradw(const Vector &x, Vector &f){
    f(0) = -16*scale*x(0)*x(1);
    f(1) = -8*scale*pow(x(0), 2);
}

//Boundary values for psi
double boundary_psi(const Vector &x){
    return scale*pow(x(0), 4)*x(1);
}

void boundary_gradpsi(const Vector &x, Vector &f){
    f(0) = 4*scale*pow(x(0), 3)*x(1);
    f(1) = scale*pow(x(0), 4);
}

//Scaling for the boundary conditions
double left_border(const Vector &x){
    double width = out_rad - int_rad;
    return 0.5*(1 - tanh((5*border/width)*(int_rad + width/border - x(0))));
}

double right_border(const Vector &x){
    double width = out_rad - int_rad;
    return 0.5*(1 - tanh((5*border/width)*(x(0) - out_rad + width/border)));
}

double lower_border(const Vector &x){
    return 0.5*(1 - tanh((5*border/height)*(height/border - x(1))));
}

double upper_border(const Vector &x){
    return 0.5*(1 - tanh((5*border/height)*(x(1) - height + height/border)));
}

double left_grad(const Vector &x){
    double width = out_rad - int_rad;
    return (2.5*border/width)*(1 - pow(tanh((5*border/width)*(int_rad + width/border - x(0))), 2));
}

double right_grad(const Vector &x){
    double width = out_rad - int_rad;
    return -(2.5*border/width)*(1 - pow(tanh((5*border/width)*(x(0) - out_rad + width/border)), 2));
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
