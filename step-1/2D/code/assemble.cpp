#include "header.h"

//Right hand side of the equation
double rhs(const Vector &x){                 
    double r_2 = pow(x(0),2);
    double z_2 = pow(x(1),2);
    return x(0)*(pow(10,-4)/7)*(2*(pow(height,2) - z_2) + pow(out_rad,2) - r_2);
}

//for Neumann condition
double boundary(const Vector &x){            
    double z_2 = pow(x(1),2);
    return x(0)*(pow(10,-4)/7)*int_rad*(pow(height,2) - z_2);
}

//Exact solution for the equation
double exact(const Vector &x){               
    double r_2 = pow(x(0),2);
    double z_2 = pow(x(1),2);
    return (pow(10,-4)/7)*0.5*(pow(height,2) - z_2)*(pow(out_rad,2) - r_2);
}

double rf(const Vector &x){
    return x(0);
}

void r_hatf(const Vector &x, Vector &f){
    f(0) = 1;
    f(1) = 0;
}

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list;

    //Dirchlet(essential) boundary conditions
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1] = ess_bdr[3] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Neumann boundary conditions
    Array<int> nbc_marker(pmesh->bdr_attributes.Max());
    nbc_marker = 0;
    nbc_marker[2] = 1;

    //Define biliniar form
    a = new ParBilinearForm(fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator(r));
    //a->AddDomainIntegrator(new ConvectionIntegrator(r_hat));
    a->Assemble();

    //Create RHS
    b = new ParLinearForm(fespace);
    FunctionCoefficient boundarycoeff(boundary);
    b->AddBoundaryIntegrator(new BoundaryLFIntegrator(boundarycoeff), nbc_marker);
    FunctionCoefficient f(rhs);
    b->AddDomainIntegrator(new DomainLFIntegrator(f));
    b->Assemble();

    //Define solution x
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(u);

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
