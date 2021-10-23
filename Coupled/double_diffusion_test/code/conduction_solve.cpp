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
    M_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K

    //Create LHS
    delete T_theta;
    delete T_e_theta;
    T_theta = Add(1., *M_0_theta, scaled_dt, *K_0_theta);
    T_e_theta = T_theta->EliminateRowsCols(ess_tdof_theta);

    delete T_phi;
    delete T_e_phi;
    T_phi = Add(1., *M_0_phi, scaled_dt, *K_0_phi);
    T_e_phi = T_phi->EliminateRowsCols(ess_tdof_phi);

    delete T;
    T = new BlockOperator(block_true_offsets);
    T->SetDiagonalBlock(0, T_theta);
    T->SetDiagonalBlock(1, T_phi);

    delete T_e;
    T_e = new BlockOperator(block_true_offsets);
    T_e->SetDiagonalBlock(0, T_e_theta);
    T_e->SetDiagonalBlock(1, T_e_phi);

    T_theta->GetDiag(T_d.GetBlock(0));
    T_phi->GetDiag(T_d.GetBlock(1));
    
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

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
    T_solver.Mult(Z, X_new);

    return 0;
}
