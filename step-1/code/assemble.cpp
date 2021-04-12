#include "header.h"

double rhs(const Vector &x){
    double r_2 = pow(x(0),2) + pow(x(1),2);
    double z_2 = pow(x(2),2);
    return (pow(10,-4)/7)*(2*(pow(height,2) - z_2) + pow(out_rad,2) - r_2);
}

double boundary(const Vector &x){
    double z_2 = pow(x(2),2);
    return (pow(10,-4)/7)*int_rad*(pow(height,2) - z_2);
}

double exact(const Vector &x){
    double r_2 = pow(x(0),2) + pow(x(1),2);
    double z_2 = pow(x(2),2);
    return (pow(10,-4)/7)*0.5*(z_2 - pow(height,2))*(r_2 - pow(out_rad,2));
}

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list;

    //Dirchlet(essential) boundary conditions
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[0] = ess_bdr[2] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Neumann boundary conditions
    Array<int> nbc_marker(pmesh->bdr_attributes.Max());
    nbc_marker = 0;
    nbc_marker[3] = 1;

    //Define biliniar form
    a = new ParBilinearForm(fespace);
    ConstantCoefficient one(1.);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
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
    u = new FunctionCoefficient(exact);
    x->ProjectCoefficient(*u);

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
