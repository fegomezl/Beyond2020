#include "header.h"
double g(const Vector &x){
    double xi(x(0));
    double yi(x(1));
    double zi(x(2));

    return 30*zi*zi*(30-zi);
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
    b->AddDomainIntegrator(new DomainLFIntegrator(one));
    FunctionCoefficient gcoeff(g);
    b->AddBoundaryIntegrator(new BoundaryLFIntegrator(gcoeff), nbc_marker);
    b->Assemble();

    //Define solution x
    x = new ParGridFunction(fespace);
    *x = 0.;

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
