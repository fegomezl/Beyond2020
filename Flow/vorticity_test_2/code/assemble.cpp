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
    GridFunctionCoefficient eta(&theta);
    ProductCoefficient neg_eta(-1., eta);

    //RHS coefficients
    FunctionCoefficient f_coeff(f_rhs);

    //Dirichlet coefficients
    FunctionCoefficient w_coeff(boundary_w);
    FunctionCoefficient psi_coeff(boundary_psi);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient neg_eta_r_inv_hat(neg_eta, r_inv_hat);
    ProductCoefficient r_f(r, f_coeff);

    //Define essential boundary conditions
    //   
    //                  1
    //            /------------\
    // (psi,w=0)  |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0

    Array<int> ess_tdof_w;
    Array<int> ess_bdr_w(pmesh->bdr_attributes.Max());
    ess_bdr_w[0] = 0; ess_bdr_w[1] = 0;
    ess_bdr_w[2] = 1; ess_bdr_w[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_w);

    Array<int> ess_tdof_psi;
    Array<int> ess_bdr_psi(pmesh->bdr_attributes.Max());
    ess_bdr_psi[0] = 1; ess_bdr_psi[1] = 1;
    ess_bdr_psi[2] = 1; ess_bdr_psi[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_psi);

    //Define grid functions
    HypreParVector W(fespace);
    w =  new ParGridFunction(fespace);
    w->ProjectCoefficient(w_coeff);
    w->ParallelProject(W);

    HypreParVector Psi(fespace);
    psi = new ParGridFunction(fespace);
    psi->ProjectCoefficient(psi_coeff);
    psi->ParallelProject(Psi);

    //Define the RHS
    HypreParVector B_w(fespace);
    ParLinearForm g(fespace);
    g.Assemble();
    g.ParallelAssemble(B_w);

    HypreParVector B_psi(fespace);
    ParLinearForm f(fespace);
    f.AddDomainIntegrator(new DomainLFIntegrator(r_f));
    f.Assemble();
    f.ParallelAssemble(B_psi);

    //Define bilinear forms of the system
    ParBilinearForm m(fespace);
    m.AddDomainIntegrator(new MassIntegrator);
    m.Assemble();
    //m.EliminateEssentialBC(ess_bdr_w, *w, g, Operator::DIAG_ONE);
    m.Finalize();
    M = m.ParallelAssemble();

    ParBilinearForm d(fespace);
    d.AddDomainIntegrator(new DiffusionIntegrator(neg_eta));
    d.AddDomainIntegrator(new ConvectionIntegrator(neg_eta_r_inv_hat));
    d.Assemble();
    //d.EliminateEssentialBC(ess_bdr_psi, *psi, f, Operator::DIAG_KEEP);
    d.Finalize();
    D = d.ParallelAssemble();

    ParMixedBilinearForm c(fespace, fespace);
    c.AddDomainIntegrator(new MixedGradGradIntegrator);
    c.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c.Assemble();
    //c.EliminateTrialDofs(ess_bdr_psi, *psi, g);
    //c.EliminateTestDofs(ess_bdr_w);
    c.Finalize();
    C = c.ParallelAssemble();

    ParMixedBilinearForm ct(fespace, fespace);
    ct.AddDomainIntegrator(new MixedGradGradIntegrator);
    ct.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    ct.Assemble();
    //ct.EliminateTrialDofs(ess_bdr_w, *w, f);
    //ct.EliminateTestDofs(ess_bdr_psi);
    ct.Finalize();
    Ct = ct.ParallelAssemble();

    //Eliminate essential DOFs
    HypreParMatrix *M_e = M->EliminateRowsCols(ess_tdof_w);
    C->EliminateRows(ess_tdof_w);
    HypreParMatrix *C_e = C->EliminateCols(ess_tdof_psi);

    HypreParMatrix *D_e = D->EliminateRowsCols(ess_tdof_psi);
    Ct->EliminateRows(ess_tdof_psi);
    HypreParMatrix *Ct_e = Ct->EliminateCols(ess_tdof_w);

    C_e->Mult(Psi, B_w, -1., 1.);
    EliminateBC(*M, *M_e, ess_tdof_w, W, B_w);

    Ct_e->Mult(W, B_psi, -1., 1.);
    EliminateBC(*D, *D_e, ess_tdof_psi, Psi, B_psi);

    //Set global system    
    X.GetBlock(0) = W;      X.GetBlock(1) = Psi;
    B.GetBlock(0) = B_w;    B.GetBlock(1) = B_psi;

    //Delete used memory
    delete M_e;
    delete D_e;
    delete C_e;
    delete Ct_e;
}

double r_f(const Vector &x){
    return x(0);
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0), -1);
    f(1) = 0.;
}

//Temperature field
double temperature_f(const Vector &x){
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = (Zmax + Zmin)/2;
    double sigma = 1;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (abs(x(0) - mid_x) > sigma && abs(x(1) - mid_y) < sigma)
        return -10;
    else
        return 10;
}

//Right hand side of the equation
double f_rhs(const Vector &x){                 
    return 0;
}

//Boundary values for w
double boundary_w(const Vector &x){
    return 0.;
}

//Boundary values for psi
double boundary_psi(const Vector &x){
    return -0.5*pow(x(0), 2);
}
