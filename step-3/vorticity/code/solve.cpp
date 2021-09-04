#include "header.h"

void rot_f(const Vector &x, DenseMatrix &f);

void Artic_sea::solve_system(){

    //Create the complete bilinear operator:
    //
    //   H = [ M     C ] 
    //       [ C^t   D ] 
    Array2D<HypreParMatrix*> hBlocks(2,2);
    hBlocks = NULL;
    hBlocks(0, 0) = M;
    hBlocks(0, 1) = C;
    hBlocks(1, 0) = Ct;
    hBlocks(1, 1) = D;

    Array2D<double> blockCoeff(2,2);
    blockCoeff(0, 0) = 1.;
    blockCoeff(0, 1) = 1.;
    blockCoeff(1, 0) = 1.;
    blockCoeff(1, 1) = 1.;

    HypreParMatrix *H = HypreParMatrixFromBlocks(hBlocks, &blockCoeff);

    SuperLUSolver *superlu = new SuperLUSolver(MPI_COMM_WORLD);
    Operator *SLU_A = new SuperLURowLocMatrix(*H);
    superlu->SetOperator(*SLU_A);
    superlu->SetPrintStatistics(true);
    superlu->SetSymmetricPattern(true);
    superlu->SetColumnPermutation(superlu::PARMETIS);
    superlu->SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    superlu->Mult(B, X);

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
    delete superlu;
    delete SLU_A;
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = -1.;
    f(1,0) = 1.;  f(1,1) = 0.;
}
