#include "header.h"

double rhs(const Vector &x){
    double r = sqrt(pow(x(0),2)+pow(x(1),2));
    double z = x(2);
    double f = (2 - out_rad/r)*(3*height - 2*z)*pow(z,2) +
               3*(r - 2*out_rad)*(height - 2*z)*r;
    return f;
}

double exact(const Vector &x){
    double r = sqrt(pow(x(0),2)+pow(x(1),2));
    double z = x(2);
    double f = 0.5*(2*out_rad - r)*r*(3*height - 2*z)*pow(z,2);
    return f;
}

double g(const Vector &x){
    double xi(x(0));
    double yi(x(1));
    double zi(x(2));

    return 30*zi*zi*(30-zi);
}

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list(0);
    Array<int> nbc_marker;
    if (pmesh->bdr_attributes.Size()){
        //dirchlet(essential) boundary conditions
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 0;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        //Neumann boundary conditions
        nbc_marker.SetSize(pmesh->bdr_attributes.Max());
        nbc_marker = 0;
        nbc_marker[12] = nbc_marker[11] = 1;
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
