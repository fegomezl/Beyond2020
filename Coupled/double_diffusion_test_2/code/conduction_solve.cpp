#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt

    //Initialize the corresponding vectors
    Z = 0.;
    dX_dt = 0.;

    //Set up RHS
    K_0->Mult(X, Z);
    for (int ii = 0; ii < ess_tdof_list.Size(); ii++)
        Z(ess_tdof_list[ii]) = 0.;

    //Solve the system
    M_solver->Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K

    //Create LHS
    delete T_e;
    HypreParMatrix *T_aux;
    T_aux = Add(1., *M_0, scaled_dt, *K_0);
    T_e = T_aux->EliminateRowsCols(ess_tdof_list);
    T_aux->AssembleDiagonal(T_d);

    delete T;
    T = new SuperLURowLocMatrix(*T_aux);
    delete T_aux;

    //Set solver
    //if (T_solver)
    //    T_solver->DismantleGrid();
    delete T_solver;
    T_solver = new SuperLUSolver(MPI_COMM_WORLD);
    T_solver->SetPrintStatistics(false);
    T_solver->SetSymmetricPattern(true);
    T_solver->SetColumnPermutation(superlu::PARMETIS);
    T_solver->SetIterativeRefine(superlu::SLU_DOUBLE);
    T_solver->SetOperator(*T);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){    
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    //Initialize the corresponding vectors
    Z = 0.;
    X_new = X;

    //Set up RHS
    M_0->Mult(X, Z);
    T_e->Mult(X, dZ);
    Z -= dZ;
    double jj = 0;
    for (int ii = 0; ii < ess_tdof_list.Size(); ii++){
        jj = ess_tdof_list[ii];
        Z(jj) = X(jj)*T_d(jj);
    } 

    //Solve the system
    T_solver->Mult(Z, X_new);

    return 0;
}
