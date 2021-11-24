#include "header.h"

void Flow_Operator::Solve(BlockVector &Z, Vector &V, Vector &rV){
    //Create block operator
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
    BlockOperator H(block_true_offsets);
    H.SetBlock(0,0,M);
    H.SetBlock(0,1,C);
    H.SetBlock(1,0,Ct);
    H.SetBlock(1,1,D);

    solver.SetOperator(H);
    prec.SetOperator(H);

    //Create the complete RHS
    BlockVector B(block_true_offsets);
    B.GetBlock(0) = B_w;
    B.GetBlock(1) = B_psi;

    //Solve the linear system Ax=B
    solver.Mult(B, Z);

    //Recover the solution on each proccesor
    w.Distribute(Z.GetBlock(0));
    psi.Distribute(Z.GetBlock(1));

    //Calculate velocity
    ParGridFunction psi_grad(&fespace_v);
    grad.Mult(psi, psi_grad);

    ParGridFunction rv(&fespace_v);
    VectorGridFunctionCoefficient Psi_grad(&psi_grad);
    MatrixVectorProductCoefficient rot_Psi_grad(rot, Psi_grad);
    rv.ProjectDiscCoefficient(rot_Psi_grad, GridFunction::ARITHMETIC);
    rv.ParallelProject(rV);

    ParGridFunction v(&fespace_v);
    FunctionCoefficient inv_R(inv_r);
    ScalarVectorProductCoefficient coeff_V(inv_R, rot_Psi_grad);
    v.ProjectDiscCoefficient(coeff_V, GridFunction::ARITHMETIC);
    v.ParallelProject(V);
}
