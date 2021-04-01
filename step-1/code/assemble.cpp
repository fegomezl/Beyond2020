#include "header.h"

void Artic_sea::assemble_system(){
    //Set the boundary values to Dirichlet
    Array<int> ess_tdof_list;
    if (pmesh->bdr_attributes.Size()){
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    //Define biliniar form
    a = new ParBilinearForm(fespace);
    ConstantCoefficient one(1.);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();

    //Create RHS
    b = new ParLinearForm(fespace);
    b->AddDomainIntegrator(new DomainLFIntegrator(one));
    b->Assemble();

    //Define solution x
    x = new ParGridFunction(fespace);
    *x = 0.;

    //Create the linear system Ax=B
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);
}
