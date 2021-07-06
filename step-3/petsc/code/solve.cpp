#include "header.h"
//#include "petsc.h"

void Artic_sea::solve_system(){
  
  PetscParMatrix *_A = new PetscParMatrix(pmesh->GetComm(), A, Operator::PETSC_MATAIJ);
  PetscLinearSolver *solver = new PetscLinearSolver(pmesh->GetComm());
  //PetscPreconditioner *prec = new PetscPreconditioner(*_A,"solver");
  PetscPreconditioner *prec = new PetscFieldSplitSolver(pmesh->GetComm(),*_A,"prec_");
  solver->SetOperator(*_A);
  //solver->SetPreconditioner(*prec);
  solver->SetAbsTol(1.e-6);
  solver->SetRelTol(1.e-6);
  solver->SetMaxIter(5000);
  solver->SetPrintLevel(2);
  solver->Mult(B, X);

  /*mfem::PetscParMatrix *pmat= new PetscParMatrix(pmesh->GetComm(),A, mfem::Operator::PETSC_MATAIJ);
  //set the local block size of the matrix
  Mat sub;
  // pmat->operator  Mat();
  //MatNestGetSubMat(*pmat,0,0,&sub);
  // MatSetBlockSize(sub,dim);
  mfem::PetscPreconditioner *prec = new mfem::PetscFieldSplitSolver(pmesh->GetComm(),*pmat,"prec_");
  mfem::PetscLinearSolver *psol = new mfem::PetscLinearSolver(pmesh->GetComm());
  psol->SetOperator(*pmat);
  psol->SetPreconditioner(*prec);
  psol->SetAbsTol(1.e-6);
  psol->SetRelTol(1.e-6);
  psol->SetMaxIter(5000);
  psol->SetPrintLevel(2);
  psol->Mult(B, X);
  //psol->GetConverged();*/

     //Recover the solution on each proccesor
    w->MakeRef(fespace_w, x.GetBlock(0), 0);
    w->Distribute(&(X.GetBlock(0)));

    psi->MakeRef(fespace_psi, x.GetBlock(1), 0);
    psi->Distribute(&(X.GetBlock(1)));

    v = new ParGridFunction(fespace_v);
    GradientGridFunctionCoefficient psi_grad(psi);
    v->ProjectCoefficient(psi_grad);

    //Delete used memory
    delete _A, solver, prec;
}
