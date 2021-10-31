#include "header.h"

void rot_f(const Vector &x, DenseMatrix &f);

double inv_r(const Vector &x);

void Artic_sea::solve_system(){

    //Create block operator
    BlockOperator H(block_true_offsets);
    H.SetBlock(0,0,M);
    H.SetBlock(0,1,C);
    H.SetBlock(1,0,Ct);
    H.SetBlock(1,1,D);

    // Solver using PETSc
    PetscLinearSolver solver(MPI_COMM_WORLD);
    PetscFieldSplitSolver prec(MPI_COMM_WORLD, H,"prec_");
    solver.SetOperator(H);
    solver.SetPreconditioner(prec);
    solver.SetAbsTol(0.0);
    solver.SetTol(1e-12);
    solver.SetMaxIter(300);
    solver.SetPrintLevel(0);
    solver.Mult(B, X);

    /*//Create block preconditioner
    HypreParVector M_d(MPI_COMM_WORLD, M->GetGlobalNumRows(), M->GetRowStarts());
    M->GetDiag(M_d);
    HypreBoomerAMG M_solver(*M);
    M_solver.iterative_mode = false;
    M_solver.SetPrintLevel(0);

    HypreParMatrix S(*C);    
    S.InvScaleRows(M_d);
    S = *ParMult(Ct, &S);
    S = *Add(-1., *D, 1., S);
    HypreBoomerAMG S_solver(S);
    S_solver.iterative_mode = false;
    S_solver.SetPrintLevel(0);

    //BlockDiagonalPreconditioner P(block_true_offsets);
    BlockLowerTriangularPreconditioner P(block_true_offsets);
    P.SetBlock(1,0,Ct);
    P.SetDiagonalBlock(0, &M_solver);
    P.SetDiagonalBlock(1, &S_solver);

    //Solve system
    GMRESSolver solver(MPI_COMM_WORLD);
    solver.SetAbsTol(0);
    solver.SetRelTol(1E-8);
    solver.SetMaxIter(1000);
    solver.SetOperator(H);
    solver.SetPreconditioner(P);
    solver.SetPrintLevel(1);
    solver.Mult(B, X);*/

    //Recover the solution on each proccesor
    w->Distribute(X.GetBlock(0));    
    psi->Distribute(X.GetBlock(1));

    //Calculate velocity
    v = new ParGridFunction(fespace_v);
    ParGridFunction psi_grad(fespace_v);

    DiscreteLinearOperator grad(fespace, fespace_v);
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();
    grad.Mult(*psi, psi_grad);

    MatrixFunctionCoefficient rot(dim, rot_f);
    FunctionCoefficient Inv_r(inv_r);

    VectorGridFunctionCoefficient Psi_grad(&psi_grad);
    MatrixVectorProductCoefficient rV(rot, Psi_grad);
    ScalarVectorProductCoefficient V(Inv_r, rV);
    v->ProjectDiscCoefficient(V, GridFunction::ARITHMETIC);
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = -1.;
    f(1,0) = 1.;  f(1,1) = 0.;
}

double inv_r(const Vector &x){
    return pow(x(0), -1);
}
