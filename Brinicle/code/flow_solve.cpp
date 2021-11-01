#include "header.h"

void Flow_Operator::Solve(BlockVector &Z, Vector &V, Vector &rV){
    //Create the complete bilinear operator:
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
    Array2D<HypreParMatrix*> HBlocks(2,2);
    HBlocks(0, 0) = M;
    HBlocks(0, 1) = C;
    HBlocks(1, 0) = Ct;
    HBlocks(1, 1) = D;

    HypreParMatrix *H = HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(*H);

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
    w.Distribute(Y.GetBlock(0));
    psi.Distribute(Y.GetBlock(1));

    //Calculate velocity
    grad.Mult(psi, psi_grad);
    Psi_grad.SetGridFunction(&psi_grad);
    rot_Psi_grad.SetBCoef(Psi_grad);
    FunctionCoefficient inv_R(inv_r);
    ScalarVectorProductCoefficient coeff_V(inv_R, rot_Psi_grad);
    v.ProjectDiscCoefficient(coeff_V, GridFunction::ARITHMETIC);
    rv.ProjectDiscCoefficient(rot_Psi_grad, GridFunction::ARITHMETIC);

    //Export flow information
    w.ParallelAverage(Z.GetBlock(0));
    psi.ParallelAverage(Z.GetBlock(1));
    v.ParallelAverage(V);
    rv.ParallelAverage(rV);

    delete H;
}
