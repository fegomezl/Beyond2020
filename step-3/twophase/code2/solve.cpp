#include "header.h"

void Artic_sea::solve_system(){
    //Create preconditioner objects
    HypreParVector *Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
                                            M->GetRowStarts());
    HypreParMatrix *Md_inv_Ct = NULL;
    HypreParMatrix *S = NULL;

    M->GetDiag(*Md);
    Md_inv_Ct = C->Transpose();
    Md_inv_Ct->InvScaleRows(*Md);
    S = ParMult(C, Md_inv_Ct);              //S = C diag(M)^-1 C^t

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
    cout.setstate(ios_base::failbit);
    solver.Mult(B, X);
    cout.clear();

    //Recover the solution on each proccesor
    v = new ParGridFunction;
    v->MakeRef(fespace_rt, x.GetBlock(0), 0);
    v->Distribute(&(X.GetBlock(0)));

    p = new ParGridFunction;
    p->MakeRef(fespace_l2, x.GetBlock(1), 0);
    p->Distribute(&(X.GetBlock(1)));

    //Prepare rules for the error analysis
    int order_quad = max(2, 2*config.order + 1);
    const IntegrationRule *irs[Geometry::NumGeom];
    for (int ii = 0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    //Calculate errors
    v_error = (v->ComputeL2Error(v_exact_coeff))/ComputeGlobalLpNorm(2, v_exact_coeff, *pmesh, irs);
    p_error = (p->ComputeL2Error(p_exact_coeff))/ComputeGlobalLpNorm(2, p_exact_coeff, *pmesh, irs);

    //Delete used memory
    delete Md, Md_inv_Ct, S;
    delete M_inv, S_inv;
    delete Prec;
}
