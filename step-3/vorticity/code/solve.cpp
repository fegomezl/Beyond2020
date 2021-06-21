#include "header.h"

void Artic_sea::solve_system(){
    //Set the preconditioner for w
    HypreBoomerAMG *amg_w = new HypreBoomerAMG(A_w);
    amg_w->SetPrintLevel(0);

    //Solve the linear system Ax=B for w
    HyprePCG *pcg_w = new HyprePCG(A_w);
    pcg_w->SetPreconditioner(*amg_w);
    pcg_w->SetPrintLevel(0);
    pcg_w->SetTol(1e-12);
    pcg_w->SetMaxIter(200);
    *X_w = 0.;
    pcg_w->Mult(*B_w, *X_w);

    //Recover the solution on each proccesor
    a_w->RecoverFEMSolution(*X_w, *b_w, *x_w);
    
    //Read the solution w for the rhs on psi
    GridFunctionCoefficient w_coeff(x_w);
    b_psi->AddDomainIntegrator(new DomainLFIntegrator(w_coeff));
    b_psi->Assemble();
    a_psi->FormLinearSystem(ess_tdof_list_psi, *x_psi, *b_psi, A_psi, *X_psi, *B_psi);

    //Set the preconditioner for psi
    HypreBoomerAMG *amg_psi = new HypreBoomerAMG(A_psi);
    amg_psi->SetPrintLevel(0);

    //Solve the linear system Ax=B for psi
    HyprePCG *pcg_psi =  new HyprePCG(A_psi);
    pcg_psi->SetPreconditioner(*amg_psi);
    pcg_psi->SetPrintLevel(0);
    pcg_psi->SetTol(1e-12);
    pcg_psi->SetMaxIter(200);
    *X_psi = 0.;
    pcg_psi->Mult(*B_psi, *X_psi);

    //Recover the solution on each proccesor
    a_psi->RecoverFEMSolution(*X_psi, *b_psi, *x_psi);

    //Recover the gradient of psi
    GradientGridFunctionCoefficient grad_psi(x_psi);
    x_v = new ParGridFunction(fespace_v);
    x_v->ProjectCoefficient(grad_psi);

    //Delete used memory
    delete amg_w, amg_psi;
    delete pcg_w, amg_psi;
}
