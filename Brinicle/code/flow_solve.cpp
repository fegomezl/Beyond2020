#include "header.h"

void Flow_Operator::Solve(BlockVector &Z, Vector &V, Vector &rV){
    //Create the complete bilinear operator
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
    Array2D<HypreParMatrix*> HBlocks(2,2);
    HBlocks(0, 0) = M;
    HBlocks(0, 1) = C;
    HBlocks(1, 0) = Ct;
    HBlocks(1, 1) = D;

    HypreParMatrix H = *HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(H);

    //Create the complete RHS
    BlockVector B(block_true_offsets);
    B.GetBlock(0) = B_w;
    B.GetBlock(1) = B_psi;

    //Set solver
    SuperLUSolver superlu(MPI_COMM_WORLD);
    superlu.SetOperator(A);
    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    superlu.Mult(B, Z);
    superlu.DismantleGrid();

    //Recover the solution on each proccesor
    w.Distribute(Z.GetBlock(0));
    psi.Distribute(Z.GetBlock(1));

    //Calculate velocity
    grad.Mult(psi, psi_grad);
    Psi_grad.SetGridFunction(&psi_grad);
    rot_Psi_grad.SetBCoef(Psi_grad);
    rv.ProjectDiscCoefficient(rot_Psi_grad, GridFunction::ARITHMETIC);
    rv.ParallelProject(rV);

    FunctionCoefficient inv_R(inv_r);
    ScalarVectorProductCoefficient coeff_V(inv_R, rot_Psi_grad);
    v.ProjectDiscCoefficient(coeff_V, GridFunction::ARITHMETIC);
    v.ParallelProject(V);
}
