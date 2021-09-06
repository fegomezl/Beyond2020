#include "header.h"

void Flow_Operator::Solve(Vector &W, Vector &Psi, Vector &V){
    //Create the complete bilinear operator:
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
    hBlocks(0, 0) = M;
    hBlocks(0, 1) = C;
    hBlocks(1, 0) = Ct;
    hBlocks(1, 1) = D;

    if (H) delete H;
    H = HypreParMatrixFromBlocks(hBlocks, &blockCoeff);

    if (SLU_A) delete SLU_A;
    SLU_A = new SuperLURowLocMatrix (*H);
    superlu.SetOperator(*SLU_A);

    //Solve the linear system Ax=B
    superlu.Mult(B, X);
    superlu.DismantleGrid();

    //Recover the solution on each proccesor
    w.Distribute(&(X.GetBlock(0)));
    psi.Distribute(&(X.GetBlock(1)));

    //Calculate velocity
    grad.Mult(psi, psi_grad);
    Psi_grad.SetGridFunction(&psi_grad);
    rot_Psi_grad.SetBCoef(Psi_grad);
    v.ProjectDiscCoefficient(rot_Psi_grad, GridFunction::ARITHMETIC);

    //Export flow information
    w.ParallelAverage(W);
    psi.ParallelAverage(Psi);
    v.ParallelAverage(V);
}
