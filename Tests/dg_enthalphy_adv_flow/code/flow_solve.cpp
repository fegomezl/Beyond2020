#include "header.h"

void Flow_Operator::Solve(BlockVector &Y, Vector &Velocity, Vector &r_Velocity){
    //Create the complete bilinear operator
    //
    //   H = [ A00  A01 ]
    //       [ A10  A11 ]
    Array2D<HypreParMatrix*> HBlocks(2,2);
    HBlocks(0, 0) = A00;
    HBlocks(0, 1) = A01;
    HBlocks(1, 0) = A10;
    HBlocks(1, 1) = A11;

    HypreParMatrix H = *HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(H);

    //Create the complete RHS
    BlockVector B(block_offsets_H1);
    B.GetBlock(0) = B0;
    B.GetBlock(1) = B1;

    //Set solver
    SuperLUSolver superlu(MPI_COMM_WORLD);
    superlu.SetOperator(A);
    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    superlu.Mult(B, Y);
    superlu.DismantleGrid();

    //Recover the solution on each proccesor
    vorticity.Distribute(Y.GetBlock(0));
    stream.Distribute(Y.GetBlock(1));

    //Calculate velocity
    ParGridFunction stream_gradient(&fespace_ND);
    gradient.Mult(stream, stream_gradient);

    ParGridFunction r_velocity(&fespace_ND);
    VectorGridFunctionCoefficient coeff_stream_gradient(&stream_gradient);
    MatrixVectorProductCoefficient coeff_rot_stream_gradient(coeff_rot, coeff_stream_gradient);
    r_velocity.ProjectDiscCoefficient(coeff_rot_stream_gradient, GridFunction::ARITHMETIC);
    r_velocity.ParallelProject(r_Velocity);

    ParGridFunction velocity(&fespace_ND);
    ScalarVectorProductCoefficient coeff_velocity(coeff_r_inv, coeff_rot_stream_gradient);
    velocity.ProjectDiscCoefficient(coeff_velocity, GridFunction::ARITHMETIC);
    velocity.ParallelProject(Velocity);
}
