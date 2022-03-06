#include "header.h"

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt

    //Initialize the corresponding vectors
    HypreParVector dTheta_dt(&fespace), dPhi_dt(&fespace);
    HypreParVector Theta(&fespace),     Phi(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        Theta(ii - block_true_offsets[0]) = X(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        Phi(ii - block_true_offsets[1]) = X(ii);
    Z_theta = 0.;   Z_phi = 0.;
    dTheta_dt = 0.; dPhi_dt = 0.; 
    dX_dt = 0.;

    //Set up RHS
    K_0_theta->Mult(-1., Theta, 1., Z_theta);   
    Z_theta += F_theta;   
    EliminateBC(*M_theta, *M_e_theta, ess_tdof_theta, dTheta_dt, Z_theta);

    K_0_phi->Mult(-1., Phi, 1., Z_phi);
    Z_phi += F_phi;
    EliminateBC(*M_phi, *M_e_phi, ess_tdof_phi, dPhi_dt, Z_phi);

    //Solve the system
    M_theta_solver.Mult(Z_theta, dTheta_dt); M_phi_solver.Mult(Z_phi, dPhi_dt);

    //Recover solution on block vector
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        dX_dt(ii) = dTheta_dt(ii - block_true_offsets[0]);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        dX_dt(ii) = dPhi_dt(ii - block_true_offsets[1]);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K

    //Create LHS
    if (T_theta) delete T_theta;
    if (T_e_theta) delete T_e_theta;
    T_theta = Add(1., *M_0_theta, scaled_dt, *K_0_theta);
    T_e_theta = T_theta->EliminateRowsCols(ess_tdof_theta);
    T_theta_prec.SetOperator(*T_theta);
    T_theta_solver.SetOperator(*T_theta);

    if (T_phi) delete T_phi;
    if (T_e_phi) delete T_e_phi;
    T_phi = Add(1., *M_0_phi, scaled_dt, *K_0_phi);
    T_e_phi = T_phi->EliminateRowsCols(ess_tdof_phi);
    T_phi_prec.SetOperator(*T_phi);
    T_phi_solver.SetOperator(*T_phi);

    //Create RHS
    dt_F_theta.Set(scaled_dt, F_theta); dt_F_phi.Set(scaled_dt, F_phi);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){    
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    //Initialize the corresponding vectors
    HypreParVector Theta_new(&fespace), Phi_new(&fespace);
    HypreParVector Theta(&fespace),     Phi(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        Theta(ii - block_true_offsets[0]) = X(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        Phi(ii - block_true_offsets[1]) = X(ii);
    Z_theta = 0.;   Z_phi = 0.;
    Theta_new = Theta; Phi_new = Phi; 
    X_new = X;

    //Set up RHS
    M_0_theta->Mult(Theta, Z_theta);   
    Z_theta += dt_F_theta;   
    EliminateBC(*T_theta, *T_e_theta, ess_tdof_theta, Theta_new, Z_theta);

    M_0_phi->Mult(Phi, Z_phi);
    Z_phi += dt_F_phi;
    EliminateBC(*T_phi, *T_e_phi, ess_tdof_phi, Phi_new, Z_phi);

    //Solve the system
    T_theta_solver.Mult(Z_theta, Theta_new); T_phi_solver.Mult(Z_phi, Phi_new);

    //Recover solution on block vector
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        X_new(ii) = Theta_new(ii - block_true_offsets[0]);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        X_new(ii) = Phi_new(ii - block_true_offsets[1]);

    return 0;
}
