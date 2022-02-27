#include "header.h"

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = B
    //Solve M(dX_dt) + K(X) = B for dX_dt

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "a\n";
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Initialize the corresponding vectors
    HypreParVector dX0_dt(&fespace_L2), dX1_dt(&fespace_L2);
    HypreParVector X0(&fespace_L2),     X1(&fespace_L2);
    for (int ii = block_offsets_L2[0]; ii < block_offsets_L2[1]; ii++)
        X0(ii - block_offsets_L2[0]) = X(ii);
    for (int ii = block_offsets_L2[1]; ii < block_offsets_L2[2]; ii++)
        X1(ii - block_offsets_L2[1]) = X(ii);
    Z0 = 0.;   Z1 = 0.;
    dX0_dt = 0.; dX1_dt = 0.;
    dX_dt = 0.;

    //Set up RHS
    K0->Mult(-1., X0, 1., Z0);
    Z0.Add(1., *B0);

    K1->Mult(-1., X1, 1., Z1);
    Z1.Add(1., *B1);

    //Solve the system  
    M0_solver.Mult(Z0, dX0_dt); M1_solver.Mult(Z1, dX1_dt); 

    //Recover solution on block vector          
    for (int ii = block_offsets_L2[0]; ii < block_offsets_L2[1]; ii++)
        dX_dt(ii) = dX0_dt(ii - block_offsets_L2[0]);   
    for (int ii = block_offsets_L2[1]; ii < block_offsets_L2[2]; ii++)
        dX_dt(ii) = dX1_dt(ii - block_offsets_L2[1]);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K
    
    if (T0) delete T0;
    T0 = Add(1., *M0, scaled_dt, *K0);
    T0_prec.SetOperator(*T0);
    T0_solver.SetOperator(*T0);

    if (T1) delete T1;
    T1 = Add(1., *M1, scaled_dt, *K1);
    T1_prec.SetOperator(*T1);
    T1_solver.SetOperator(*T1);

    //Set dt for RHS
    if (B0_dt) delete B0_dt;
    B0_dt = new HypreParVector(&fespace_L2);
    B0_dt->Set(scaled_dt, *B0);

    if (B1_dt) delete B1_dt;
    B1_dt = new HypreParVector(&fespace_L2);
    B1_dt->Set(scaled_dt, *B1);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = B
    //Solve M(X_new - X) + dt*K(X_new) = dt*B for X_new
    
    //Initialize the corresponding vectors
    HypreParVector X0_new(&fespace_L2), X1_new(&fespace_L2);
    HypreParVector X0(&fespace_L2),     X1(&fespace_L2);
    for (int ii = block_offsets_L2[0]; ii < block_offsets_L2[1]; ii++)
        X0(ii - block_offsets_L2[0]) = X(ii);
    for (int ii = block_offsets_L2[1]; ii < block_offsets_L2[2]; ii++)
        X1(ii - block_offsets_L2[1]) = X(ii);
    Z0 = 0.;   Z1 = 0.;
    X0_new = X0; X1_new = X1;
    X_new = X;

    //Set up RHS
    M0->Mult(X0, Z0);
    Z0.Add(1., *B0_dt);

    M1->Mult(X1, Z1);
    Z1.Add(1., *B1_dt);

    //Solve the system  
    T0_solver.Mult(Z0, X0_new); T1_solver.Mult(Z1, X1_new); 
    //Recover solution on block vector          
    for (int ii = block_offsets_L2[0]; ii < block_offsets_L2[1]; ii++)
        X_new(ii) = X0_new(ii - block_offsets_L2[0]);   
    for (int ii = block_offsets_L2[1]; ii < block_offsets_L2[2]; ii++)
        X_new(ii) = X1_new(ii - block_offsets_L2[1]);

    return 0;
}
