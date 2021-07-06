#include "header.h"

void Artic_sea::solve_system(){

    Array2D<HypreParMatrix*> hBlocks(2,2);
    hBlocks = NULL;
    hBlocks(0, 0) = M;
    hBlocks(0, 1) = C;
    hBlocks(1, 0) = Ct;
    hBlocks(1, 1) = D;

    Array2D<double> blockCoeff(2,2);
    blockCoeff = 1.0;

    HypreParMatrix *H = HypreParMatrixFromBlocks(hBlocks, &blockCoeff);

    SuperLUSolver *superlu = new SuperLUSolver(MPI_COMM_WORLD);
    Operator *SLU_A = new SuperLURowLocMatrix(*H);
    superlu->SetOperator(*SLU_A);
    superlu->SetPrintStatistics(true);
    superlu->SetSymmetricPattern(true);

    int slu_colperm = 4;
    int slu_rowperm = 1;
    int slu_iterref = 2;

    if (slu_colperm == 0)
      {
	superlu->SetColumnPermutation(superlu::NATURAL);
      }
    else if (slu_colperm == 1)
      {
	superlu->SetColumnPermutation(superlu::MMD_ATA);
      }
    else if (slu_colperm == 2)
      {
	superlu->SetColumnPermutation(superlu::MMD_AT_PLUS_A);
      }
    else if (slu_colperm == 3)
      {
	superlu->SetColumnPermutation(superlu::COLAMD);
      }
    else if (slu_colperm == 4)
      {
	superlu->SetColumnPermutation(superlu::METIS_AT_PLUS_A);
      }
    else if (slu_colperm == 5)
      {
	superlu->SetColumnPermutation(superlu::PARMETIS);
      }
    else if (slu_colperm == 6)
      {
	superlu->SetColumnPermutation(superlu::ZOLTAN);
      }
    
    if (slu_rowperm == 0)
      {
	superlu->SetRowPermutation(superlu::NOROWPERM);
      }
    else if (slu_rowperm == 1)
      {
#ifdef MFEM_USE_SUPERLU5
	superlu->SetRowPermutation(superlu::LargeDiag);
#else
	superlu->SetRowPermutation(superlu::LargeDiag_MC64);
#endif
      }
    
    if (slu_iterref == 0)
      {
	superlu->SetIterativeRefine(superlu::NOREFINE);
      }
    else if (slu_iterref == 1)
      {
	superlu->SetIterativeRefine(superlu::SLU_SINGLE);
      }
    else if (slu_iterref == 2)
      {
	superlu->SetIterativeRefine(superlu::SLU_DOUBLE);
      }
    else if (slu_iterref == 3)
      {
	superlu->SetIterativeRefine(superlu::SLU_EXTRA);
      }

    //Solve the linear system Ax=B
    X = 0.;
    superlu->Mult(B, X);

    //Recover the solution on each proccesor
    w->Distribute(&(X.GetBlock(0)));
    
    psi->Distribute(&(X.GetBlock(1)));

    v = new ParGridFunction(fespace_v);
    GradientGridFunctionCoefficient psi_grad(psi);
    v->ProjectCoefficient(psi_grad);

    //Delete used memory
    delete H;
    delete superlu;
    delete SLU_A;
}
