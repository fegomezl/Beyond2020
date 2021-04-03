#include "header.h"

double rhs(const Vector &x){
    double r_2 = pow(x(0),2)+pow(x(1),2);
    double z_2 = pow(x(2),2);
    return 2*(z_2 - pow(height,2)) + r_2 - pow(out_rad,2);
}

double g(const Vector &x){
    double z_2 = pow(x(2),2);
    return int_rad*(pow(height,2) - z_2);
}

double exact(const Vector &x){
    double r_2 = pow(x(0),2)+pow(x(1),2);
    double z_2 = pow(x(2),2);
    return 0.5*(pow(height,2) - z_2)*(r_2 - pow(out_rad,2));
}

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list;
    Array<int> nbc_marker;
    if (pmesh->bdr_attributes.Size()){
        //dirchlet(essential) boundary conditions
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        //Neumann boundary conditions
        nbc_marker.SetSize(pmesh->bdr_attributes.Max());
        nbc_marker = 0;
        nbc_marker[12] = 1;
    }

    //Define biliniar form
    a = new ParBilinearForm(fespace);
    ConstantCoefficient one(1.);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();

    //Create RHS
    b = new ParLinearForm(fespace);
    FunctionCoefficient gcoeff(g);
    b->AddBoundaryIntegrator(new BoundaryLFIntegrator(gcoeff), nbc_marker);
    FunctionCoefficient f(rhs);
    b->AddDomainIntegrator(new DomainLFIntegrator(f));
    b->Assemble();

    //Define solution x
    x = new ParGridFunction(fespace);
    u = new FunctionCoefficient(exact);
    //x->ProjectCoefficient(*u);

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
