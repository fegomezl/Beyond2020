#include "header.h"

//From  M(dX_dt) + K(X) = B
//Solve M(dX_dt) + K(X) = B for dX_dt
void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    
    //Initialize the corresponding vectors
    HypreParVector dX0_dt(&fespace_H1), dX1_dt(&fespace_H1);
    HypreParVector X0(&fespace_H1),     X1(&fespace_H1);
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        X0(ii - block_offsets_H1[0]) = X(ii);
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        X1(ii - block_offsets_H1[1]) = X(ii);
    Z0 = 0.;   Z1 = 0.;
    dX0_dt = 0.; dX1_dt = 0.;

    //Set up RHS
    K0->Mult(-1., X0, 1., Z0);
    Z0.Add(1., B0);

    K1->Mult(-1., X1, 1., Z1);
    Z1.Add(1., B1);

    //Solve the system  
    M0_solver.Mult(Z0, dX0_dt); M1_solver.Mult(Z1, dX1_dt); 

    //Recover solution on block vector          
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        dX_dt(ii) = dX0_dt(ii - block_offsets_H1[0]);   
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        dX_dt(ii) = dX1_dt(ii - block_offsets_H1[1]);
}

// Setup the system T*X = M*RHS + scaled_dt*B
// With T = M + scaled+dt*K
int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    
    //Create matrix T = M + scaled+dt*K
    if (T0) delete T0;
    T0 = Add(1., *M0, scaled_dt, *K0);
    T0_prec.SetOperator(*T0);
    T0_solver.SetOperator(*T0);

    if (T1) delete T1;
    T1 = Add(1., *M1, scaled_dt, *K1);
    T1_prec.SetOperator(*T1);
    T1_solver.SetOperator(*T1);

    //Set dt for RHS
    B0_dt.Set(scaled_dt, B0);
    B1_dt.Set(scaled_dt, B1);

    *j_status = 1;
    return 0;
}

//From  M(dX_dt) + K(X) = B
//Solve M(X - RHS) + scaled_dt*K(X) = scaled_dt*B for X
int Transport_Operator::SUNImplicitSolve(const Vector &RHS, Vector &X, double tol){
    
    //Initialize the corresponding vectors
    HypreParVector X0(&fespace_H1),   X1(&fespace_H1);
    HypreParVector RHS0(&fespace_H1), RHS1(&fespace_H1);
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++){
        RHS0(ii - block_offsets_H1[0]) = RHS(ii);
        X0(ii - block_offsets_H1[0]) = X(ii);
    }
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++){
        RHS1(ii - block_offsets_H1[1]) = RHS(ii);
        X1(ii - block_offsets_H1[1]) = X(ii);
    }
    Z0 = 0.;   Z1 = 0.;

    //Set up RHS
    M0->Mult(X0, Z0);
    Z0.Add(1., B0_dt);

    M1->Mult(X1, Z1);
    Z1.Add(1., B1_dt);

    //Solve the system  
    T0_solver.Mult(Z0, X0); T1_solver.Mult(Z1, X1); 

    //Recover solution on block vector          
    for (int ii = block_offsets_H1[0]; ii < block_offsets_H1[1]; ii++)
        X(ii) = X0(ii - block_offsets_H1[0]);   
    for (int ii = block_offsets_H1[1]; ii < block_offsets_H1[2]; ii++)
        X(ii) = X1(ii - block_offsets_H1[1]);

    return 0;
}
