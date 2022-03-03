#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(X, *rV);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Z, *V, *rV);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        theta->Distribute(X.GetBlock(0));
        phi->Distribute(X.GetBlock(1));
        w->Distribute(Z.GetBlock(0));
        psi->Distribute(Z.GetBlock(1));
        v->Distribute(V);
        rv->Distribute(rV);

        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
        }

        //Normalize stream
        if (config.rescale){
            double psi_local_max = psi->Max(), psi_max;
            double psi_local_min = psi->Min(), psi_min;
            MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            for (int ii = 0; ii < psi->Size(); ii++)
                    (*psi)(ii) = ((*psi)(ii)-psi_min)/(psi_max-psi_min);
        }

        //Graph
        paraview_out->SetCycle(vis_impressions);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }

    //Print the system state
    double percentage = 100*(t-config.t_init)/(config.t_final-config.t_init);
    string progress = to_string((int)percentage)+"%";
    if (config.master){
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << dt << setw(12)
             << t  << setw(12)
             << progress << "\r";
        cout.flush();

        std::ofstream out;
        out.open("results/progress.txt", std::ios::app);
        out << left << setw(12)
            << iteration << setw(12)
            << dt << setw(12)
            << t  << setw(12)
            << progress << "\n";
        out.close();
    }
}

void Conduction_Operator::SetParameters(const BlockVector &X, const Vector &rV){
    //Recover actual information
    ParGridFunction theta(&fespace_dg), phi(&fespace_dg), phase(&fespace_dg), rv(&fespace_v);
    theta.Distribute(X.GetBlock(0));
    phi.Distribute(X.GetBlock(1));
    rv.Distribute(rV); 

    //Associate the values of each auxiliar function
    double DT = 0.;
    for (int ii = 0; ii < phase.Size(); ii++){
        DT = theta(ii) - T_fun(phi(ii));
        theta(ii) = DT;
        phase(ii) = 0.5*(1 + tanh(5*config.invDeltaT*DT));
    }

    ParGridFunction aux_C(phase), aux_K(phase), aux_D(phase), aux_L(phase);
    aux_C *= config.c_l-config.c_s; aux_C += config.c_s;
    aux_K *= config.k_l-config.k_s; aux_K += config.k_s;
    aux_D *= config.d_l-config.d_s; aux_D += config.d_s;
    aux_L *= config.L_l-config.L_s; aux_L += config.L_s;

    aux_C.ExchangeFaceNbrData();
    aux_K.ExchangeFaceNbrData();
    aux_D.ExchangeFaceNbrData();
    aux_L.ExchangeFaceNbrData();

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);
    GridFunctionCoefficient coeff_L(&aux_L);

    //Construct latent heat term
    GradientGridFunctionCoefficient dT(&theta), dH(&phase);
    InnerProductCoefficient dHdT(dH, dT), dT_2(dT, dT);

    SumCoefficient dT_2e(config.EpsilonT, dT_2);
    
    PowerCoefficient inv_dT_2(dT_2e, -1.);
    ProductCoefficient DeltaT(dHdT, inv_dT_2);
    ProductCoefficient LDeltaT(coeff_L, DeltaT);

    SumCoefficient coeff_CL(coeff_C, LDeltaT);

    //Construct final coefficients
    ProductCoefficient coeff_rC(coeff_r, coeff_CL);
    ProductCoefficient coeff_rK(coeff_r, coeff_K);
    ProductCoefficient coeff_rD(coeff_r, coeff_D);

    VectorGridFunctionCoefficient coeff_rV(&rv);
    ScalarVectorProductCoefficient coeff_rCV(coeff_CL, coeff_rV);

    //Create corresponding bilinear forms
    delete M_theta;
    ParBilinearForm m_theta(&fespace_dg);
    m_theta.AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m_theta.Assemble();
    m_theta.Finalize();
    M_theta = m_theta.ParallelAssemble();

    M_theta_prec.SetOperator(*M_theta);
    M_theta_solver.SetOperator(*M_theta);

    delete K_theta;
    ParBilinearForm k_theta(&fespace_dg);
    k_theta.AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k_theta.AddDomainIntegrator(new ConvectionIntegrator(coeff_rCV));
    k_theta.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rK, config.sigma, config.kappa));
    k_theta.AddInteriorFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rCV, 1.));
    k_theta.AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rK, config.sigma, config.kappa), ess_bdr_theta);
    k_theta.AddBdrFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rCV, 1.), ess_bdr_theta);
    k_theta.Assemble();
    k_theta.Finalize();
    K_theta = k_theta.ParallelAssemble();

    delete K_phi;
    ParBilinearForm k_phi(&fespace_dg);
    k_phi.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi.AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k_phi.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rD, config.sigma, config.kappa));
    k_phi.AddInteriorFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.));
    k_phi.AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rD, config.sigma, config.kappa), ess_bdr_phi);
    k_phi.AddBdrFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.), ess_bdr_phi);
    k_phi.Assemble();
    k_phi.Finalize();
    K_phi = k_phi.ParallelAssemble();

    if (B_theta) delete B_theta;
    ConstantCoefficient theta_nu(theta_n);
    ParLinearForm b_theta(&fespace_dg);
    b_theta.AddBdrFaceIntegrator(new DGDirichletLFIntegrator(theta_nu, coeff_rK, config.sigma, config.kappa), ess_bdr_theta);
    b_theta.AddBdrFaceIntegrator(new BoundaryFlowIntegrator(theta_nu, coeff_rCV, -1.), ess_bdr_theta);
    b_theta.Assemble();
    B_theta = b_theta.ParallelAssemble();

    if (B_phi) delete B_phi;
    ConstantCoefficient phi_nu(phi_n);
    ParLinearForm b_phi(&fespace_dg);
    b_phi.AddBdrFaceIntegrator(new DGDirichletLFIntegrator(phi_nu, coeff_rD, config.sigma, config.kappa), ess_bdr_phi);
    b_phi.AddBdrFaceIntegrator(new BoundaryFlowIntegrator(phi_nu, coeff_rV, -1.), ess_bdr_phi);
    b_phi.Assemble();
    B_phi = b_phi.ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){

    //Recover actual information
    ParGridFunction theta(&fespace_dg), phi(&fespace_dg);
    ParGridFunction theta_dr(&fespace), phi_dr(&fespace);
    ParGridFunction eta(&fespace_dg);

    theta.Distribute(X.GetBlock(0));
    phi.Distribute(X.GetBlock(1));

    theta.GetDerivative(1, 0, theta_dr);
    phi.GetDerivative(1, 0, phi_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < eta.Size(); ii++){
        double T = theta(ii);
        double S = phi(ii);
        double P = 0.5*(1 + tanh(5*config.invDeltaT*(theta(ii) - T_fun(phi(ii))))); 

        eta(ii) = config.EpsilonEta + pow(1-P, 2)/(pow(P, 3) + config.EpsilonEta);

        theta(ii) = delta_rho_t_fun(T, S);
        phi(ii) = delta_rho_p_fun(T, S);
    }

    theta.ExchangeFaceNbrData();
    phi.ExchangeFaceNbrData();
    eta.ExchangeFaceNbrData();

    //Properties coefficients
    GridFunctionCoefficient Eta(&eta);
    ProductCoefficient neg_Eta(-1., Eta);

    GridFunctionCoefficient Theta_dr(&theta_dr);
    GridFunctionCoefficient k_t(&theta);
    ProductCoefficient k_Theta_dr(k_t, Theta_dr);

    GridFunctionCoefficient Phi_dr(&phi_dr);
    GridFunctionCoefficient k_p(&phi);
    ProductCoefficient k_Phi_dr(k_p, Phi_dr);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient neg_Eta_r_inv_hat(neg_Eta, r_inv_hat);
    ProductCoefficient k_r_Theta_dr(coeff_r, k_Theta_dr);
    ProductCoefficient k_r_Phi_dr(coeff_r, k_Phi_dr);

    //Define non-constant bilinear forms of the system
    delete D;
    ParBilinearForm d(&fespace);
    d.AddDomainIntegrator(new DiffusionIntegrator(neg_Eta));
    d.AddDomainIntegrator(new ConvectionIntegrator(neg_Eta_r_inv_hat));
    d.Assemble();
    d.Finalize();
    D = d.ParallelAssemble();
    D_e = D->EliminateRowsCols(ess_tdof_psi);

    //Define the non-constant RHS
    ParLinearForm f(&fespace);
    f.AddDomainIntegrator(new DomainLFIntegrator(k_r_Theta_dr));
    f.AddDomainIntegrator(new DomainLFIntegrator(k_r_Phi_dr));
    f.Assemble();
    f.ParallelAssemble(B_psi);
    Ct_e->Mult(W, B_psi, -1., 1.);
    EliminateBC(*D, *D_e, ess_tdof_psi, Psi, B_psi);
}
