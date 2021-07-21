#include "header.h"

void rot_f(const Vector &x, DenseMatrix &f);

int Artic_sea::solve_system(){

    //Create the complete bilinear operator:
    //
    //   H = [ M    C ] 
    //       [ C^t  cm D ] 
    Array2D<HypreParMatrix*> hBlocks(2,2);
    hBlocks = NULL;
    hBlocks(0, 0) = M;
    hBlocks(0, 1) = C;
    hBlocks(1, 0) = C->Transpose();
    hBlocks(1, 1) = D;

    Array2D<double> blockCoeff(2,2);
    blockCoeff(0, 0) = 1.;
    blockCoeff(0, 1) = -1.;
    blockCoeff(1, 0) = -1.;
    blockCoeff(1, 1) = -1.;

     HypreParMatrix *H = HypreParMatrixFromBlocks(hBlocks, &blockCoeff);
    // std::unique_ptr<HypreParMatrix> H(new HypreParMatrixFromBlocks(hBlocks, &blockCoeff));

    //std::unique_ptr<SuperLUSolver> superlu = std::make_unique<SuperLUSolver>(MPI_COMM_WORLD);
    //std::unique_ptr<SuperLURowLocMatrix> SLU_A = std::make_unique<SuperLURowLocMatrix>(H);
    SuperLUSolver *superlu = new SuperLUSolver(MPI_COMM_WORLD);
    SuperLURowLocMatrix *SLU_A = new SuperLURowLocMatrix(*H);

    superlu->SetOperator(*SLU_A);
    superlu->SetPrintStatistics(true);
    superlu->SetSymmetricPattern(true);
    superlu->SetColumnPermutation(superlu::PARMETIS);
    superlu->SetIterativeRefine(superlu::SLU_DOUBLE);

    X.Randomize();
    superlu->Mult(B, X);
    superlu ->DismantleGrid();


    /*PetscParMatrix *_A = new PetscParMatrix(pmesh->GetComm(), H, Operator::PETSC_MATAIJ);
    PetscLinearSolver *solver = new PetscLinearSolver(pmesh->GetComm());
    //PetscLinearSolver* solver = new PetscLinearSolver(pmesh->GetComm());
    //PetscPreconditioner* prec = new PetscPreconditioner(A,"solver");
    PetscPreconditioner *prec = new PetscFieldSplitSolver(pmesh->GetComm(),*_A,"prec_");
    solver->SetOperator(*_A);
    //solver->SetPreconditioner(*prec);
    solver->SetAbsTol(1.e-6);
    solver->SetRelTol(1.e-6);
    solver->SetMaxIter(5000);
    solver->SetPrintLevel(1);
    solver->Mult(B, X);

    PetscLinearSolver * petsc = new PetscLinearSolver(MPI_COMM_WORLD);
    // Convert to PetscParMatrix
    petsc->SetOperator(PetscParMatrix(H, Operator::PETSC_MATAIJ));
    petsc->SetPrintLevel(2);
    petsc->Mult(B, X);*/


    //Recover the solution on each proccesor
    w->Distribute(&(X.GetBlock(0)));
    for (int ii = 0; ii < w->Size(); ii++)
        (*w)(ii) += (*w_aux)(ii); 
    
    psi->Distribute(&(X.GetBlock(1)));
    for (int ii = 0; ii < psi->Size(); ii++)
        (*psi)(ii) += (*psi_aux)(ii); 

    v = new ParGridFunction(fespace_v);
    GradientGridFunctionCoefficient psi_grad(psi);
    MatrixFunctionCoefficient rot(dim, rot_f);
    MatrixVectorProductCoefficient rV(rot, psi_grad);
    v->ProjectDiscCoefficient(rV, GridFunction::ARITHMETIC);

    //Delete used memory
    delete H;
    delete superlu;
    delete SLU_A;

    return 0;
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = 1.;
    f(1,0) = -1.; f(1,1) = 0.;
}
