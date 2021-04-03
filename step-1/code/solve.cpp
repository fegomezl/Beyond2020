#include "header.h"

void Artic_sea::solve_system(){
    //Set the preconditioner
    HypreBoomerAMG *amg = new HypreBoomerAMG(A);
    amg->SetPrintLevel(0);

    //Solve the linear system Ax=B
    HyprePCG *pcg = new HyprePCG(A);
    pcg->SetPreconditioner(*amg);
    pcg->SetPrintLevel(0);
    pcg->SetTol(1e-12);
    pcg->SetMaxIter(200);
    pcg->Mult(B, X);

    //Recover the solution on each proccesor
    a->RecoverFEMSolution(X, *b, *x);

    //Calculate error
    l2_error = x->ComputeL2Error(*u);
    double volume = M_PI*height*(pow(out_rad,2) - pow(int_rad,2));
    l2_error /= volume;

    //Delete used memory
    delete amg;
    delete pcg;
}
