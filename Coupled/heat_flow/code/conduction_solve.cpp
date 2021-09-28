#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt
    
    Z = 0.;
    dX_dt = 0.;
    
    K_0->Mult(-1., X, 1., Z);
    EliminateBC(*M, *M_e, ess_tdof_list, dX_dt, Z);

    M_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K

    if (T) delete T;
    if (T_e) delete T_e;
    T = Add(1., *M_0, scaled_dt, *K_0);
    T_e = T->EliminateRowsCols(ess_tdof_list);
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    Z = 0.;
    X_new = X;

    M_0->Mult(X, Z);
    EliminateBC(*T, *T_e, ess_tdof_list, X_new, Z);

    T_solver.Mult(Z, X_new);
    return 0;
}
