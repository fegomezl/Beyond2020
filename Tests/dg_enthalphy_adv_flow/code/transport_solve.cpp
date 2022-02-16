#include "header.h"

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = B
    //Solve M(dX_dt) + K(X) = B for dX_dt
    Z = 0.;
    dX_dt = 0.;
    
    K->Mult(-1., X, 1., Z);
    Z.Add(1., *B);

    M_solver.Mult(Z, dX_dt);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K
    if (T) delete T;
    T = Add(1., *M, scaled_dt, *K);
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

    //Set dt for RHS
    if (B_dt) delete B_dt;
    B_dt = new HypreParVector(&fespace_L2);
    B_dt->Set(scaled_dt, *B);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = B
    //Solve M(X_new - X) + dt*K(X_new) = dt*B for X_new
    Z = 0.;
    X_new = X;

    M->Mult(X, Z);
    Z.Add(1., *B_dt);

    T_solver.Mult(Z, X_new);
    return 0;
}
