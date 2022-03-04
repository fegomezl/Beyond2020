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

    HypreParMatrix H = *HypreParMatrixFromBlocks(HBlocks);
    SuperLURowLocMatrix A(H);

    //Create the complete RHS
    B.GetBlock(0) = G;
    B.GetBlock(1) = F;

    SuperLUSolver superlu(MPI_COMM_WORLD);
    superlu.SetOperator(A);
    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    Z = 0.;
    superlu.Mult(B, Z);
    superlu.DismantleGrid();

    //Calculate velocity
    psi.Distribute(Z.GetBlock(1));
    grad.Mult(psi, psi_grad);
    Psi_grad.SetGridFunction(&psi_grad);
    rot_Psi_grad.SetBCoef(Psi_grad);
    FunctionCoefficient inv_R(r_inv_f);
    ScalarVectorProductCoefficient coeff_V(inv_R, rot_Psi_grad);
    v.ProjectDiscCoefficient(coeff_V, GridFunction::ARITHMETIC);
    rv.ProjectDiscCoefficient(rot_Psi_grad, GridFunction::ARITHMETIC);

    //Export flow information
    v.ParallelAverage(V);
    rv.ParallelAverage(rV);
}
