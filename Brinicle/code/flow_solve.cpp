#include "header.h"

//Solution of the current system
void Flow_Operator::Solve(BlockVector &Y, Vector &Velocity, Vector &rVelocity){

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
    B.GetBlock(0) = *B0;
    B.GetBlock(1) = *B1;

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

    //Calculate velocity field
    gradient.Mult(stream, stream_gradient);
    VectorGridFunctionCoefficient coeff_stream_gradient(&stream_gradient);
    MatrixVectorProductCoefficient coeff_rV(coeff_rot, coeff_stream_gradient);
    ScalarVectorProductCoefficient coeff_V(coeff_r_inv, coeff_rV);
    velocity.ProjectDiscCoefficient(coeff_V, GridFunction::ARITHMETIC);
    rvelocity.ProjectDiscCoefficient(coeff_rV, GridFunction::ARITHMETIC);

    //Export flow information
    velocity.ParallelAverage(Velocity);
    rvelocity.ParallelAverage(rVelocity);

    delete H;
}
