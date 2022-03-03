#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt

    //Initialize the corresponding vectors
    HypreParVector dTheta_dt(&fespace_dg), dPhi_dt(&fespace_dg);
    HypreParVector Theta(&fespace_dg),     Phi(&fespace_dg);
    for (int ii = block_true_offsets_dg[0]; ii < block_true_offsets_dg[1]; ii++)
        Theta(ii - block_true_offsets_dg[0]) = X(ii);
    for (int ii = block_true_offsets_dg[1]; ii < block_true_offsets_dg[2]; ii++)
        Phi(ii - block_true_offsets_dg[1]) = X(ii);
    Z_theta = 0.;   Z_phi = 0.;
    dTheta_dt = 0.; dPhi_dt = 0.; 
    dX_dt = 0.;

    //Set up RHS
    K_theta->Mult(-1., Theta, 1., Z_theta);   
    Z_theta += *B_theta;

    K_phi->Mult(-1., Phi, 1., Z_phi);
    Z_phi += *B_phi;

    //Solve the system
    M_theta_solver.Mult(Z_theta, dTheta_dt); M_phi_solver.Mult(Z_phi, dPhi_dt);

    //Recover solution on block vector
    for (int ii = block_true_offsets_dg[0]; ii < block_true_offsets_dg[1]; ii++)
        dX_dt(ii) = dTheta_dt(ii - block_true_offsets_dg[0]);
    for (int ii = block_true_offsets_dg[1]; ii < block_true_offsets_dg[2]; ii++)
        dX_dt(ii) = dPhi_dt(ii - block_true_offsets_dg[1]);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K

    //Create LHS
    if (T_theta) delete T_theta;
    T_theta = Add(1., *M_theta, scaled_dt, *K_theta);
    T_theta_prec.SetOperator(*T_theta);
    T_theta_solver.SetOperator(*T_theta);

    if (T_phi) delete T_phi;
    T_phi = Add(1., *M_phi, scaled_dt, *K_phi);
    T_phi_prec.SetOperator(*T_phi);
    T_phi_solver.SetOperator(*T_phi);

    if (B_dt_theta) delete B_dt_theta;
    B_dt_theta = new HypreParVector(&fespace_dg);
    B_dt_theta->Set(scaled_dt, *B_theta);
    
    if (B_dt_phi) delete B_dt_phi;
    B_dt_phi = new HypreParVector(&fespace_dg);
    B_dt_phi->Set(scaled_dt, *B_phi);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){    
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    //Initialize the corresponding vectors
    HypreParVector Theta_new(&fespace_dg), Phi_new(&fespace_dg);
    HypreParVector Theta(&fespace_dg),     Phi(&fespace_dg);
    for (int ii = block_true_offsets_dg[0]; ii < block_true_offsets_dg[1]; ii++)
        Theta(ii - block_true_offsets_dg[0]) = X(ii);
    for (int ii = block_true_offsets_dg[1]; ii < block_true_offsets_dg[2]; ii++)
        Phi(ii - block_true_offsets_dg[1]) = X(ii);
    Z_theta = 0.;   Z_phi = 0.;
    Theta_new = Theta; Phi_new = Phi; 
    X_new = X;

    //Set up RHS
    M_theta->Mult(Theta, Z_theta);   
    Z_theta += *B_dt_theta;

    M_phi->Mult(Phi, Z_phi);
    Z_phi += *B_dt_phi;

    //Solve the system
    T_theta_solver.Mult(Z_theta, Theta_new); T_phi_solver.Mult(Z_phi, Phi_new);

    //Recover solution on block vector
    for (int ii = block_true_offsets_dg[0]; ii < block_true_offsets_dg[1]; ii++)
        X_new(ii) = Theta_new(ii - block_true_offsets_dg[0]);
    for (int ii = block_true_offsets_dg[1]; ii < block_true_offsets_dg[2]; ii++)
        X_new(ii) = Phi_new(ii - block_true_offsets_dg[1]);

    return 0;
}
