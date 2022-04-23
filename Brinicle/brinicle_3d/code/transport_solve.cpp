#include "header.h"

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = B
    //Solve M(dX_dt) + K(X) = B for dX_dt
    
    //Initialize the corresponding vectors
    HypreParVector dX0_dt(&fespace_H1), dX1_dt(&fespace_H1);
    HypreParVector X0(&fespace_H1),     X1(&fespace_H1);
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        X0(ii - block_offsets_H1[0]) = X(ii);
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        X1(ii - block_offsets_H1[1]) = X(ii);
    Z0 = 0.;   Z1 = 0.;
    dX0_dt = 0.; dX1_dt = 0.;
    dX_dt = 0.;

    //Set up RHS
    K0->Mult(-1., X0, 1., Z0);
    Z0.Add(1., *B0);
    EliminateBC(*M0, *M0_e, ess_tdof_0, dX0_dt, Z0);

    K1->Mult(-1., X1, 1., Z1);
    Z1.Add(1., *B1);
    EliminateBC(*M1, *M1_e, ess_tdof_1, dX1_dt, Z1);

    //Solve the system  
    M0_solver.Mult(Z0, dX0_dt); M1_solver.Mult(Z1, dX1_dt); 

    //Recover solution on block vector          
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        dX_dt(ii) = dX0_dt(ii - block_offsets_H1[0]);   
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        dX_dt(ii) = dX1_dt(ii - block_offsets_H1[1]);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K
    
    if (T0) delete T0;
    if (T0_e) delete T0_e;
    T0 = Add(1., *M0_o, scaled_dt, *K0);
    T0_e = T0->EliminateRowsCols(ess_tdof_0);
    T0_prec.SetOperator(*T0);
    T0_solver.SetOperator(*T0);

    if (T1) delete T1;
    if (T1_e) delete T1_e;
    T1 = Add(1., *M1_o, scaled_dt, *K1);
    T1_e = T1->EliminateRowsCols(ess_tdof_1);
    T1_prec.SetOperator(*T1);
    T1_solver.SetOperator(*T1);

    //Set dt for RHS
    B0_dt.Set(scaled_dt, *B0);

    B1_dt.Set(scaled_dt, *B1);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = B
    //Solve M(X_new - X) + dt*K(X_new) = dt*B for X_new
    
    //Initialize the corresponding vectors
    HypreParVector X0_new(&fespace_H1), X1_new(&fespace_H1);
    HypreParVector X0(&fespace_H1),     X1(&fespace_H1);
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        X0(ii - block_offsets_H1[0]) = X(ii);
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        X1(ii - block_offsets_H1[1]) = X(ii);
    Z0 = 0.;   Z1 = 0.;
    X0_new = X0; X1_new = X1;
    X_new = X;

    //Set up RHS
    M0_o->Mult(X0, Z0);
    Z0.Add(1., B0_dt);
    EliminateBC(*T0, *T0_e, ess_tdof_0, X0_new, Z0);

    M1_o->Mult(X1, Z1);
    Z1.Add(1., B1_dt);
    EliminateBC(*T1, *T1_e, ess_tdof_1, X1_new, Z1);

    //Solve the system  
    T0_solver.Mult(Z0, X0_new); T1_solver.Mult(Z1, X1_new); 

    //Recover solution on block vector          
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        X_new(ii) = X0_new(ii - block_offsets_H1[0]);   
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        X_new(ii) = X1_new(ii - block_offsets_H1[1]);

    return 0;
}
