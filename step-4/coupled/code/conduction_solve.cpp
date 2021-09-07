#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    HypreParVector Theta(&fespace);
    HypreParVector Phi(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        Theta(ii) = X(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        Phi(ii - block_true_offsets[1]) = X(ii);

    //Set all the system parts
    ParGridFunction theta(&fespace), phi(&fespace);
    ParLinearForm z_theta(&fespace), z_phi(&fespace);

    theta.SetFromTrueDofs(Theta); phi.SetFromTrueDofs(Phi);
    k_theta->Mult(theta, z_theta); k_phi->Mult(phi, z_phi);
    z_theta.Neg(); z_phi.Neg();

    OperatorHandle A_theta, A_phi;
    HypreParVector Z_theta(&fespace), Z_phi(&fespace);
    ParGridFunction dtheta_dt(&fespace), dphi_dt(&fespace);
    HypreParVector dTheta_dt(&fespace), dPhi_dt(&fespace);
    dtheta_dt = 0.; dphi_dt = 0.;

    //Solve the system
    m_theta->FormLinearSystem(ess_tdof_list_theta, dtheta_dt, z_theta, A_theta, dTheta_dt, Z_theta);
    m_phi->FormLinearSystem(ess_tdof_list_phi, dphi_dt, z_phi, A_phi, dPhi_dt, Z_phi);

    ParLinearForm rhs_theta(&fespace), rhs_phi(&fespace);
    HypreParVector Rhs_theta(&fespace), Rhs_phi(&fespace);

    rhs_theta.AddBoundaryIntegrator(new BoundaryLFIntegrator(r_newmann_theta), newmann_bdr_theta);
    rhs_theta.AddBoundaryIntegrator(new BoundaryLFIntegrator(neg_r_robin_h_ref_theta), robin_bdr_theta);
    rhs_theta.Assemble();
    rhs_theta.ParallelAssemble(Rhs_theta);
    Z_theta += Rhs_theta;

    rhs_phi.AddBoundaryIntegrator(new BoundaryLFIntegrator(r_newmann_phi), newmann_bdr_phi);
    rhs_phi.AddBoundaryIntegrator(new BoundaryLFIntegrator(neg_r_robin_h_ref_phi), robin_bdr_phi);
    rhs_phi.Assemble();
    rhs_phi.ParallelAssemble(Rhs_phi);
    Z_phi += Rhs_phi;

    M_theta_solver.Mult(Z_theta, dTheta_dt); M_phi_solver.Mult(Z_phi, dPhi_dt);

    //Recover solution on block vector
    dX_dt = 0.;
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        dX_dt(ii) = dTheta_dt(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        dX_dt(ii) = dPhi_dt(ii - block_true_offsets[1]);
}

//Setup the ODE Jacobian T = M + gamma*K
int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Set coefficients
    dt_coeff_rK.SetAConst(scaled_dt);           dt_coeff_rCLV.SetAConst(scaled_dt);
    dt_coeff_rD.SetAConst(scaled_dt);           dt_coeff_rV.SetAConst(scaled_dt);
    dt_r_newmann_theta.SetAConst(scaled_dt);    dt_r_newmann_phi.SetAConst(scaled_dt);
    neg_dt_r_robin_h_theta.SetAConst(scaled_dt); neg_dt_r_robin_h_phi.SetAConst(scaled_dt);
    neg_dt_r_robin_h_ref_theta.SetAConst(scaled_dt); neg_dt_r_robin_h_ref_phi.SetAConst(scaled_dt);

    //Create bilinear forms
    if(t_theta) delete t_theta;
    t_theta = new ParBilinearForm(&fespace);
    t_theta->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    t_theta->AddDomainIntegrator(new ConvectionIntegrator(dt_coeff_rCLV));
    t_theta->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rK));
    t_theta->AddBoundaryIntegrator(new MassIntegrator(neg_dt_r_robin_h_theta), robin_bdr_theta);
    t_theta->Assemble();
    t_theta->FormSystemMatrix(ess_tdof_list_theta, T_theta);
    T_theta_solver.SetOperator(T_theta);

    if(t_phi) delete t_phi;
    t_phi = new ParBilinearForm(&fespace);
    t_phi->AddDomainIntegrator(new MassIntegrator(coeff_r));
    t_phi->AddDomainIntegrator(new ConvectionIntegrator(dt_coeff_rV));
    t_phi->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rD));
    t_phi->AddBoundaryIntegrator(new MassIntegrator(neg_dt_r_robin_h_phi), robin_bdr_phi);
    t_phi->Assemble();
    t_phi->FormSystemMatrix(ess_tdof_list_phi, T_phi);
    T_phi_solver.SetOperator(T_phi);

    *j_status = 1;
    return 0;
}

//Solve the system Ax = z -> (M - gamma*K)x = Mb
int Conduction_Operator::SUNImplicitSolve(const Vector &B, Vector &X, double tol){
    //Initialize the corresponding vectors
    HypreParVector B_theta(&fespace), B_phi(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        B_theta(ii) = B(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        B_phi(ii - block_true_offsets[1]) = B(ii);

    //Set all the system parts
    HypreParVector Theta(&fespace), Phi(&fespace);
    HypreParVector Z_theta(&fespace), Z_phi(&fespace);

    //Solve the system
    M_theta.Mult(B_theta, Z_theta);      M_phi.Mult(B_phi, Z_phi);

    ParLinearForm rhs_theta(&fespace), rhs_phi(&fespace);
    HypreParVector Rhs_theta(&fespace), Rhs_phi(&fespace);

    rhs_theta.AddBoundaryIntegrator(new BoundaryLFIntegrator(dt_r_newmann_theta), newmann_bdr_theta);
    rhs_theta.AddBoundaryIntegrator(new BoundaryLFIntegrator(neg_dt_r_robin_h_ref_theta), robin_bdr_theta);
    rhs_theta.Assemble();
    rhs_theta.ParallelAssemble(Rhs_theta);
    Z_theta += Rhs_theta;

    rhs_phi.AddBoundaryIntegrator(new BoundaryLFIntegrator(dt_r_newmann_phi), newmann_bdr_phi);
    rhs_phi.AddBoundaryIntegrator(new BoundaryLFIntegrator(neg_dt_r_robin_h_ref_phi), robin_bdr_phi);
    rhs_phi.Assemble();
    rhs_phi.ParallelAssemble(Rhs_phi);
    Z_phi += Rhs_phi;

    T_theta_solver.Mult(Z_theta, Theta); T_phi_solver.Mult(Z_phi, Phi);

    //Recover solution on block vector
    X = 0.;
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        X(ii) = Theta(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        X(ii) = Phi(ii - block_true_offsets[1]);

    return 0;
}
