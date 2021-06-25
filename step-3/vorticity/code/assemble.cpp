#include "header.h"

double mid_x = (out_rad + int_rad)/2;
double mid_y = height/2;
double sigma = (out_rad - int_rad)/5;
double mag = 10;
double temp = 5;

//Temperature field
double T(const Vector &x){
    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);

    if (r_2 < sigma)
        return -temp;
    else
        return temp;
}

//Phase change fraction field
double frac(const Vector &x){
    double T_f = 0;
    double delta = 0.0001;
    return 0.5*(tanh(5*(T(x) - T_f)/delta) + 1);
}

//Permeability
double eta(const Vector &x){
    //double m = 0.00000001;
    //return pow(1 - frac(x), 2)/(pow(frac(x), 3) + m);

    /*double mag = 1e+6;
    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return mag;
    else
        return 0.1;*/
    return 0;
}

//Right hand side of the equation
double rhs(const Vector &x){                 
    double b_g = 2.06;

    return x(0)*(out_rad - int_rad)/2;
}

//Boundary for stream function
double boundary(const Vector &x){            
    return x(0);
}

double rf(const Vector &x){
    return x(0);
}

void Artic_sea::assemble_system(){
    //Dirchlet(essential) boundary conditions
    ess_bdr_w.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr_w = 0;
    ess_bdr_w[1] = ess_bdr_w[2] = ess_bdr_w[3] = 1;
    fespace_w->GetEssentialTrueDofs(ess_bdr_w, ess_tdof_list_w);

    ess_bdr_psi.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr_psi = 0;
    ess_bdr_psi[1] = ess_bdr_psi[2] = ess_bdr_psi[3] = 1;
    fespace_psi->GetEssentialTrueDofs(ess_bdr_psi, ess_tdof_list_psi);

    //Define coefficients
    ConstantCoefficient viscosity(1);
    FunctionCoefficient permeability(eta);
    FunctionCoefficient g(rhs);
    ConstantCoefficient zero(0.);
    FunctionCoefficient boundary_coeff(boundary);

    //Define biliniar forms
    a_w = new ParBilinearForm(fespace_w);
    a_w->AddDomainIntegrator(new DiffusionIntegrator(viscosity));
    a_w->AddDomainIntegrator(new MassIntegrator(permeability));
    a_w->Assemble();

    a_psi = new ParBilinearForm(fespace_psi);
    a_psi->AddDomainIntegrator(new DiffusionIntegrator());
    a_psi->Assemble();

    //Create RHS
    b_w = new ParLinearForm(fespace_w);
    b_w->AddDomainIntegrator(new DomainLFIntegrator(g));
    b_w->Assemble();

    b_psi = new ParLinearForm(fespace_psi);

    //Define solution x and apply boundary values
    x_w = new ParGridFunction(fespace_w);
    x_w->ProjectBdrCoefficient(zero, ess_bdr_w);

    x_psi = new ParGridFunction(fespace_psi);
    x_psi->ProjectBdrCoefficient(boundary_coeff, ess_bdr_psi);

    //Create True Dofs Vectors
    X_w = new HypreParVector(fespace_w);
    X_psi = new HypreParVector(fespace_psi);

    B_w = new HypreParVector(fespace_w);
    B_psi = new HypreParVector(fespace_psi);

    //Create the linear system Ax=B
    a_w->FormLinearSystem(ess_tdof_list_w, *x_w, *b_w, A_w, *X_w, *B_w);
}
