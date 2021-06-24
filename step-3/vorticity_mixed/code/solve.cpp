#include "header.h"

void Artic_sea::solve_system(){
    //Create preconditioner objects
    HypreParVector *Dd = new HypreParVector(MPI_COMM_WORLD, D->GetGlobalNumRows(),
                                            D->GetRowStarts());
    HypreParVector *Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
                                            M->GetRowStarts());
    HypreParMatrix *Md_inv_Ct = NULL;
    HypreParMatrix *C_Md_inv_Ct = NULL;
    HypreParMatrix *S = NULL;

    M->GetDiag(*Md);
    D->GetDiag(*Dd);
    Md_inv_Ct = C->Transpose();
    Md_inv_Ct->InvScaleRows(*Md);
    C_Md_inv_Ct = ParMult(C, Md_inv_Ct);              //S = C diag(M)^-1 C^t
    (*M) *= -1.;
    S = ParAdd(M, C_Md_inv_Ct);

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
    solver.SetPrintLevel(0);
    solver.SetAbsTol(1e-10);
    solver.SetRelTol(1e-6);
    solver.SetMaxIter(500);
    X = 0.;
    solver.Mult(B, X);

    //Recover the solution on each proccesor
    psi = new ParGridFunction;
    psi->MakeRef(fespace_psi, x.GetBlock(0), 0);
    psi->Distribute(&(X.GetBlock(0)));

    w = new ParGridFunction;
    w->MakeRef(fespace_w, x.GetBlock(1), 0);
    w->Distribute(&(X.GetBlock(1)));

    v = new ParGridFunction(fespace_v);
    GradientGridFunctionCoefficient psi_grad(psi);
    v->ProjectCoefficient(psi_grad);

    //Prepare rules for the error analysis
    int order_quad = max(2, 2*config.order + 1);
    const IntegrationRule *irs[Geometry::NumGeom];
    for (int ii = 0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    //Delete used memory
    delete Md, Md_inv_Ct, C_Md_inv_Ct, S;
    delete M_inv, S_inv;
    delete Prec;
}
