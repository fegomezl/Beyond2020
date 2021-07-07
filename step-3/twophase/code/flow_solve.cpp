#include "header.h"

void Flow_Operator::Solve(Config config, HypreParVector *X_Psi, ParGridFunction *x_psi, const HypreParVector *X_T){

    this->Update_T(config, X_T);

    //Create preconditioner objects
    HypreParVector *Dd = new HypreParVector(MPI_COMM_WORLD, D->GetGlobalNumRows(),
                                            D->GetRowStarts());
    //HypreParVector *Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
    //                                        M->GetRowStarts());
    HypreParMatrix *Dd_inv_Ct = NULL;
    HypreParMatrix *C_Dd_inv_Ct = NULL;
    HypreParMatrix *S = NULL;
    HypreParMatrix *Mplus = NULL;

    //M->GetDiag(*Md);
    D->GetDiag(*Dd);
    //Md_inv_Ct = C->Transpose();
    Dd_inv_Ct = Ct;
    Dd_inv_Ct->InvScaleRows(*Dd);
    C_Dd_inv_Ct = ParMult(C, Dd_inv_Ct);              //S = +/-(M - C diag(D)^-1 C^t)
    Mplus = M;
    (*Mplus) *= -1.;
    S = ParAdd(Mplus, C_Dd_inv_Ct);

    //Create solvers for the preconditioner
    Solver *D_inv = new HypreDiagScale(*D);
    D_inv->iterative_mode = false;

    Solver *S_inv = new HypreBoomerAMG(*S);
    S_inv->iterative_mode = false;

    //Create the complete preconditioner
    //
    //   P = [ diag(D)  0]
    //       [   0      S]
    BlockDiagonalPreconditioner *Prec = new BlockDiagonalPreconditioner(block_true_offsets);
    Prec->SetDiagonalBlock(0, D_inv);
    Prec->SetDiagonalBlock(1, S_inv);

    //Solve the linear system Ax=B
    MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetPreconditioner(*Prec);
    solver.SetOperator(*A);
    solver.SetPrintLevel(0);
    solver.SetAbsTol(1e-10);
    solver.SetRelTol(1e-6);
    solver.SetMaxIter(1000);
    Y = 0.;
    solver.Mult(B, Y);

    //Recover the solution on each proccesor
    psi->MakeRef(&fespace, y.GetBlock(0), 0);
    psi->Distribute(&(Y.GetBlock(0)));

    w->MakeRef(&fespace, y.GetBlock(1), 0);
    w->Distribute(&(Y.GetBlock(1)));

    //Delete used memory
    delete Dd, Mplus, Dd_inv_Ct, C_Dd_inv_Ct, S;
    delete D_inv, S_inv;
    delete Prec;
}
