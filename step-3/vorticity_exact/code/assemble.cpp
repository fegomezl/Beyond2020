#include "header.h"

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
    GridFunctionCoefficient eta(theta);

    //Define local coefficients
    ConstantCoefficient mu(config.viscosity);
    ProductCoefficient neg_mu(-1., mu);
    ProductCoefficient neg_eta(-1., eta);

    FunctionCoefficient f_coeff(f_rhs);
    ProductCoefficient neg_f_coeff(-1., f_coeff);

    FunctionCoefficient w_coeff(scaled_boundary_w);
    VectorFunctionCoefficient w_grad(dim, scaled_boundary_gradw);
    ProductCoefficient neg_mu_w(neg_mu, w_coeff);
    ScalarVectorProductCoefficient mu_w_grad(mu, w_grad);
    ScalarVectorProductCoefficient neg_mu_w_grad(neg_mu, w_grad);

    FunctionCoefficient psi_coeff(scaled_boundary_psi);
    VectorFunctionCoefficient psi_grad(dim, scaled_boundary_gradpsi);
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
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    //Define grid functions
    w =  new ParGridFunction;
    w_aux = new ParGridFunction(fespace);
    w->MakeRef(fespace, x.GetBlock(0), 0);
    w_aux->ProjectCoefficient(w_coeff);

    psi = new ParGridFunction;
    psi_aux = new ParGridFunction(fespace);
    psi->MakeRef(fespace, x.GetBlock(1), 0);
    psi_aux->ProjectCoefficient(psi_coeff);

    //Define the RHS
    g = new ParLinearForm;
    g->Update(fespace, b.GetBlock(0), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(neg_mu_w));
    g->AddDomainIntegrator(new DomainLFGradIntegrator(mu_psi_grad));
    g->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_psi_grad));
    g->Assemble();
    g->ParallelAssemble(B.GetBlock(0));

    f = new ParLinearForm;
    f->Update(fespace, b.GetBlock(1), 0);
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_f_coeff));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(mu_w_grad));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_mu_w_grad));
    f->Assemble();
    f->ParallelAssemble(B.GetBlock(1));

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace);
    m->AddDomainIntegrator(new MassIntegrator(mu));
    m->Assemble();
    m->EliminateEssentialBCFromDofs(ess_tdof_list_w, *w, *g);
    m->Finalize();
    M = m->ParallelAssemble();

    d = new ParBilinearForm(fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta));
    d->Assemble();
    d->EliminateEssentialBCFromDofs(ess_tdof_list_psi, *psi, *f);
    d->Finalize();
    D = d->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace, fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator(mu));
    c->Assemble();
    OperatorHandle Ch;
    c->FormRectangularSystemMatrix(ess_tdof_list_psi, ess_tdof_list_w, Ch);
    C = Ch.Is<HypreParMatrix>();
}

//Temperature field
double temperature_f(const Vector &x){
    double mid_x = (out_rad + int_rad)/2;
    double mid_y = height/2;
    double sigma = (out_rad - int_rad)/10;

    return 10;
    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;
    
}

//Right hand side of the equation
double f_rhs(const Vector &x){                 
  return -12*x(0)*x(0)+12*x(1)*x(1);
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 12*x(1)*x(1)-12*x(0)*x(0);
}

void boundary_gradw(const Vector &x, Vector &f){
  f(0) = -24*x(0);
  f(1) = 24*x(1);
}

//Boundary values for psi
double boundary_psi(const Vector &x){
  return pow(x(0),4)-pow(x(1),4);
}

void boundary_gradpsi(const Vector &x, Vector &f){
  f(0) = 4*pow(x(0),3);
  f(1) = -4*pow(x(1),3);
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
