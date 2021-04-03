#include "header.h"
double g(const Vector &x){
    double xi(x(0));
    double yi(x(1));
    double zi(x(2));

    return 30*zi*zi*(30-zi);
}

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

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list;
    if (pmesh->bdr_attributes.Size()){
        //Dirichlet(esscential) boundary conditions
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        ess_bdr[12] = 0;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    //Neumann boundary conditions
    Array<int> nbc_marker(pmesh->bdr_attributes.Max());
    nbc_marker = 0;
    nbc_marker[12] = 1;

    //Define biliniar form
    a = new ParBilinearForm(fespace);
    ConstantCoefficient one(1.);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();

    //Create RHS
    b = new ParLinearForm(fespace);
<<<<<<< HEAD
    b->AddDomainIntegrator(new DomainLFIntegrator(one));
    FunctionCoefficient gcoeff(g);
    b->AddBoundaryIntegrator(new BoundaryLFIntegrator(gcoeff), nbc_marker);
=======
    FunctionCoefficient f(rhs);
    b->AddDomainIntegrator(new DomainLFIntegrator(f));
>>>>>>> f4c7f9434a7405927babf912e0268e10b27c4a3c
    b->Assemble();

    //Define solution x
    x = new ParGridFunction(fespace);
    u = new FunctionCoefficient(exact);
    x->ProjectCoefficient(*u);

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
