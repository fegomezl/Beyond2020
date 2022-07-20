#include "header.h"

//Solution of the current system
void Flow_Operator::Solve(BlockVector &Y){

    //Create the complete bilinear operator:
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
    Array2D<HypreParMatrix*> HBlocks(2,2);
    HBlocks(0, 0) = A00;
    HBlocks(0, 1) = A01;
    HBlocks(1, 0) = A10;
    HBlocks(1, 1) = A11;

    HypreParMatrix *H = HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(*H);

    //Create the complete RHS
    B.GetBlock(0) = B0;
    B.GetBlock(1) = B1;

    //Create the solver object
    SuperLUSolver superlu(MPI_COMM_WORLD);
    superlu.SetOperator(A);
    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    superlu.Mult(B, Y);
    superlu.DismantleGrid();

    vorticity.Distribute(Y.GetBlock(0)); 
    stream.Distribute(Y.GetBlock(1)); 

    delete H;
}
