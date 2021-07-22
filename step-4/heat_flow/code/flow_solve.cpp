#include "header.h"

void Flow_Operator::Solve(const HypreParVector *Theta){

    this->Update_T(Theta);

    //Create the complete bilinear operator:
    //
    //   H = [ M    C ]
    //       [ C^t  D ]
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

    SuperLUSolver superlu = SuperLUSolver(MPI_COMM_WORLD);
    SuperLURowLocMatrix SLU_A(*H);
    superlu.SetOperator(SLU_A);
    superlu.SetPrintStatistics(false);
    superlu.SetSymmetricPattern(true);
    superlu.SetColumnPermutation(superlu::PARMETIS);
    superlu.SetIterativeRefine(superlu::SLU_DOUBLE);

    //Solve the linear system Ax=B
    Y.Randomize();
    superlu.Mult(B, Y);

    //Recover the solution on each proccesor
    w->Distribute(&(Y.GetBlock(0)));
    for (int ii = 0; ii < w->Size(); ii++)
        (*w)(ii) += (*w_aux)(ii);

    psi->Distribute(&(Y.GetBlock(1)));
    for (int ii = 0; ii < psi->Size(); ii++)
        (*psi)(ii) += (*psi_aux)(ii);

    //Calculate rV
    v->Randomize();
    v_aux->Randomize();

    //gradpsi.SetGridFunction(psi);
    grad.Mult(*psi, *v_aux);
    rV_aux.SetGridFunction(v_aux);
    rV.SetBCoef(rV_aux);
    v->ProjectDiscCoefficient(rV, GridFunction::ARITHMETIC);

    //Free memory
    delete H;
}
