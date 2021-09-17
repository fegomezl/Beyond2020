#include "header.h"

void rot_f(const Vector &x, DenseMatrix &f);

void Artic_sea::solve_system(){

    //Create the complete bilinear operator:
    //
    //   H = [ M     C ] 
    //       [ C^t   D ] 
    Array2D<HypreParMatrix*> HBlocks(2,2);
    HBlocks = NULL;
    HBlocks(0, 0) = M;
    HBlocks(0, 1) = C;
    HBlocks(1, 0) = Ct;
    HBlocks(1, 1) = D;

    HypreParMatrix *H = HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(*H);

    SuperLUSolver superlu(MPI_COMM_WORLD);
    superlu.SetOperator(A);
    superlu.SetPrintStatistics(true);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    superlu.Mult(B, X);

    //Recover the solution on each proccesor
    w->Distribute(&(X.GetBlock(0)));    
    psi->Distribute(&(X.GetBlock(1)));

    //Calculate velocity
    v = new ParGridFunction(fespace_v);

    ParGridFunction psi_grad(fespace_v);

    DiscreteLinearOperator grad(fespace, fespace_v);
    grad.AddDomainIntegrator(new GradientInterpolator);
    grad.Assemble();
    grad.Finalize();

    grad.Mult(*psi, psi_grad);

    VectorGridFunctionCoefficient Psi_grad(&psi_grad);
    MatrixFunctionCoefficient rot(dim, rot_f);
    MatrixVectorProductCoefficient V(rot, Psi_grad);
    v->ProjectDiscCoefficient(V, GridFunction::ARITHMETIC);

    //Delete used memory
    delete H;
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = -1.;
    f(1,0) = 1.;  f(1,1) = 0.;
}
