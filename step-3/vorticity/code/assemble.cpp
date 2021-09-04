#include "header.h"

//Rotational functions
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
    ParGridFunction theta(fespace);
    FunctionCoefficient temperature(temperature_f);
    theta.ProjectCoefficient(temperature);
    for (int ii = 0; ii < theta.Size(); ii++){
        theta(ii) = 0.5*(1 + tanh(5*config.invDeltaT*(theta(ii) - config.T_f)));
        theta(ii) = config.epsilon_eta + (1 - pow(theta(ii), 2))/(pow(theta(ii), 3) + config.epsilon_eta);
    }

    //Rotational coefficients
    FunctionCoefficient r(r_f);
    VectorFunctionCoefficient r_inv_hat(dim, r_inv_hat_f);

    //Properties coefficients
    ConstantCoefficient inv_mu(pow(config.viscosity, -1));
    GridFunctionCoefficient eta(&theta);
    ProductCoefficient neg_eta(-1., eta);

    //RHS coefficients
    FunctionCoefficient f_coeff(f_rhs);
    ProductCoefficient inv_mu_f(inv_mu, f_coeff);

    //Dirichlet coefficients
    FunctionCoefficient w_coeff(boundary_w);
    FunctionCoefficient psi_coeff(boundary_psi);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient neg_eta_r_inv_hat(neg_eta, r_inv_hat);
    ProductCoefficient inv_mu_rf(r, inv_mu_f);

    //Define essential boundary conditions
    //   
    //                  1
    //            /------------\
    // (psi,w=0)  |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0

    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w[0] = 1; ess_bdr_w[1] = 1;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 0;

    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 0;

    //Define grid functions
    w =  new ParGridFunction(fespace);
    w->ProjectCoefficient(w_coeff);

    psi = new ParGridFunction(fespace);
    psi->ProjectCoefficient(psi_coeff);

    //Define the RHS
    ParLinearForm g(fespace);
    g.Assemble();

    ParLinearForm f(fespace);
    f.AddDomainIntegrator(new DomainLFIntegrator(inv_mu_rf));
    f.Assemble();

    //Define bilinear forms of the system
    ParBilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator);
    m.Assemble();
    m.EliminateEssentialBC(ess_bdr_w, *w, g);
    m.Finalize();
    M = m.ParallelAssemble();

    ParBilinearForm d(fespace);
    d.AddDomainIntegrator(new DiffusionIntegrator(neg_eta));
    d.AddDomainIntegrator(new ConvectionIntegrator(neg_eta_r_inv_hat));
    d.Assemble();
    d.EliminateEssentialBC(ess_bdr_psi, *psi, f);
    d.Finalize();
    D = d.ParallelAssemble();

    ParMixedBilinearForm c(fespace, fespace);
    c.AddDomainIntegrator(new MixedGradGradIntegrator);
    c.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c.Assemble();

    //Forma uno
    c.EliminateTrialDofs(ess_bdr_psi, *psi, g);
    c.EliminateTestDofs(ess_bdr_w);
    c.Finalize();
    C = c.ParallelAssemble();

    ParMixedBilinearForm ct(fespace, fespace);
    ct.AddDomainIntegrator(new MixedGradGradIntegrator);
    ct.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    ct.Assemble();
    ct.EliminateTrialDofs(ess_bdr_w, *w, f);
    ct.EliminateTestDofs(ess_bdr_psi);
    ct.Finalize();
    Ct = ct.ParallelAssemble();

    //Transfer to TrueDofs
    w->ParallelAssemble(X.GetBlock(0));
    psi->ParallelAssemble(X.GetBlock(1));

    g.ParallelAssemble(B.GetBlock(0));
    f.ParallelAssemble(B.GetBlock(1));
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

    //double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return -10;
    else
        return 10;
}

//Right hand side of the equation
double f_rhs(const Vector &x){                 
    if (x(0) > 5)
        return 1;
    else
        return 0;
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 0.;
}

//Boundary values for psi
double boundary_psi(const Vector &x){
    double a = 1;
    double v = 10;
    if (x(0) < a)
        return -v*0.5*pow(x(0), 2)*x(1)/height;
    else
        return -v*0.5*pow(a, 2)*x(1)/height;
}
