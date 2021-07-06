#include "header.h"

void Artic_sea::solve_system(){
    //Create preconditioner objects
    HypreParVector *Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
                                            M->GetRowStarts());
    HypreParMatrix *Md_inv_C = NULL;
    HypreParMatrix *Ct_Md_inv_C = NULL;
    HypreParMatrix *S = NULL;

    M->GetDiag(*Md);
    Md_inv_C = C;
    Md_inv_C->InvScaleRows(*Md);
    Ct_Md_inv_C = ParMult(Ct, Md_inv_C);              //S = +/-(D - C^t diag(M)^-1 C)
    (*D) *= -1.;
    S = ParAdd(D, Ct_Md_inv_C);

    //Create solvers for the preconditioner
    Solver *M_inv = new HypreDiagScale(*M);
    M_inv->iterative_mode = false;

    Solver *S_inv = new HypreBoomerAMG(*S);
    S_inv->iterative_mode = false;

    //Create the complete preconditioner
    //
    //   P = [ diag(M)  0]
    //       [   0      S] 
    BlockDiagonalPreconditioner *Prec = new BlockDiagonalPreconditioner(block_true_offsets);
    Prec->SetDiagonalBlock(0, M_inv);
    Prec->SetDiagonalBlock(1, S_inv);

    //Solve the linear system Ax=B
    MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetPreconditioner(*Prec);
    solver.SetOperator(*A);
    solver.SetPrintLevel(1);
    solver.SetAbsTol(1e-2);
    solver.SetRelTol(1e-6);
    solver.SetMaxIter(2000);
    X = 0.;
    solver.Mult(B, X);

    //Recover the solution on each proccesor
    w->MakeRef(fespace_w, x.GetBlock(0), 0);
    w->Distribute(&(X.GetBlock(0)));
    for (int ii = 0; ii < w->Size(); ii++)
        (*w)(ii) += (*w_boundary)(ii); 

    psi->MakeRef(fespace_psi, x.GetBlock(1), 0);
    psi->Distribute(&(X.GetBlock(1)));
    for (int ii = 0; ii < psi->Size(); ii++)
        (*psi)(ii) += (*psi_boundary)(ii); 

    v = new ParGridFunction(fespace_v);
    GradientGridFunctionCoefficient psi_grad(psi);
    v->ProjectCoefficient(psi_grad);

    //Delete used memory
    delete Md, Md_inv_C, Ct_Md_inv_C, S;
    delete M_inv, S_inv;
    delete Prec;
}
