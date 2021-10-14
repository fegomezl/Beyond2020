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
        theta->Distribute(&(X.GetBlock(0)));
        phi->Distribute(&(X.GetBlock(1)));
        w->Distribute(&(Z.GetBlock(0)));
        psi->Distribute(&(Z.GetBlock(1)));
        v->Distribute(V);
        rv->Distribute(rV);
        
        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = config.T_f + T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
        }

        //Normalize stream
        double psi_local_max = psi->Max(), psi_max;
        double psi_local_min = psi->Min(), psi_min;
        MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        for (int ii = 0; ii < psi->Size(); ii++)
                (*psi)(ii) = ((*psi)(ii)-psi_min)/(psi_max-psi_min);

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
    theta.SetFromTrueDofs(X.GetBlock(0));
    phi.SetFromTrueDofs(X.GetBlock(1));
    rv.SetFromTrueDofs(rV); 

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < theta.Size(); ii++){
        double T = theta(ii) - config.T_f - T_fun(phi(ii));

        if (T > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
            aux_D(ii) = config.D_l;
        } else {
            aux_C(ii) = config.c_s;// + delta_c_s_fun(theta(ii), phi(ii));
            aux_K(ii) = config.k_s;// + delta_k_s_fun(theta(ii), phi(ii));
            aux_D(ii) = config.D_s;
        }
        aux_L(ii) = config.L;// + delta_l_s_fun(theta(ii), phi(ii));

        theta(ii) = T;
        phase(ii) = 0.5*(1 + tanh(5*config.invDeltaT*T));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);
    GridFunctionCoefficient coeff_L(&aux_L);

    //Construct latent heat term
    GradientGridFunctionCoefficient dT(&theta);
    GradientGridFunctionCoefficient dH(&phase);
    
    dHdT.SetACoef(dH);  dT_2.SetACoef(dT);
    dHdT.SetBCoef(dT);  dT_2.SetBCoef(dT);

    SumCoefficient dT_2e(config.EpsilonT, dT_2);
    
    PowerCoefficient inv_dT_2(dT_2e, -1.);
    ProductCoefficient DeltaT(dHdT, inv_dT_2);
    ProductCoefficient LDeltaT(coeff_L, DeltaT);

    SumCoefficient coeff_CL(coeff_C, LDeltaT);

    //Construct final coefficients
    coeff_rC.SetBCoef(coeff_CL);
    coeff_rK.SetBCoef(coeff_K); 
    coeff_rD.SetBCoef(coeff_D); 

    coeff_rV.SetGridFunction(&rv);
    coeff_rCV.SetACoef(coeff_CL);
    coeff_rCV.SetBCoef(coeff_rV);

    //Create corresponding bilinear forms
    if (m_theta) delete m_theta;
    if (M_theta) delete M_theta;
    if (M_e_theta) delete M_e_theta;
    if (M_0_theta) delete M_0_theta;
    m_theta = new ParBilinearForm(&fespace);
    m_theta->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m_theta->Assemble();
    m_theta->Finalize();
    M_theta = m_theta->ParallelAssemble();
    M_e_theta = M_theta->EliminateRowsCols(ess_tdof_theta);
    M_0_theta = m_theta->ParallelAssemble();

    M_theta_prec.SetOperator(*M_theta);
    M_theta_solver.SetOperator(*M_theta);

    if (m_phi) delete m_phi;
    if (M_phi) delete M_phi;
    if (M_e_phi) delete M_e_phi;
    if (M_0_phi) delete M_0_phi;
    m_phi = new ParBilinearForm(&fespace);
    m_phi->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m_phi->Assemble();
    m_phi->Finalize();
    M_phi = m_phi->ParallelAssemble();
    M_e_phi = M_phi->EliminateRowsCols(ess_tdof_phi);
    M_0_phi = m_phi->ParallelAssemble();

    M_phi_prec.SetOperator(*M_phi);
    M_phi_solver.SetOperator(*M_phi);

    if(k_theta) delete k_theta;
    if(K_0_theta) delete K_0_theta;
    k_theta = new ParBilinearForm(&fespace);
    k_theta->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k_theta->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCV));
    k_theta->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_theta), robin_bdr_theta);
    k_theta->Assemble();
    k_theta->Finalize();
    K_0_theta = k_theta->ParallelAssemble();

    if(k_phi) delete k_phi;
    if(K_0_phi) delete K_0_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k_phi->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_phi), robin_bdr_phi);
    k_phi->Assemble();
    k_phi->Finalize();
    K_0_phi = k_phi->ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Update information
    theta.SetFromTrueDofs(X.GetBlock(0));
    phi.SetFromTrueDofs(X.GetBlock(1));
    eta.SetFromTrueDofs(X.GetBlock(0));

    theta.GetDerivative(1, 0, theta_dr);
    phi.GetDerivative(1, 0, phi_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < eta.Size(); ii++){
        double T = theta(ii);
        double S = phi(ii);
        double T_f = config.T_f + T_fun(S);

        eta(ii) = 0.5*(1 + tanh(5*config.invDeltaT*(T - T_f)));
        eta(ii) = config.EpsilonEta + pow(1-eta(ii), 2)/(pow(eta(ii), 3) + config.EpsilonEta);

        double a00 = -10.27542, a01 = 0.83195,
               a10 = 0.38667,   a11 = -0.02801,
               a20 = -0.00624,  
               a30 = 0.00006;

        theta(ii) = (a00 + a10*T + a20*pow(T, 2) + a30*pow(T, 3))*S +               // k_t = g*b_t/mu
                    (a01 + a11*T)*pow(abs(S), 1.5);                                  

        double b00 = 2061.94264, b01 = -68.54083, b02 = 24.43573,
               b10 = -10.27542,  b11 = 1.24793,
               b20 = 0.19333,    b21 = -0.02101,
               b30 = -0.00208,
               b40 = 0.00001;

        phi(ii) = (b00 + b10*T + b20*pow(T, 2) + b30*pow(T, 3) + b40*pow(T, 4)) +   // k_s = g*b_s/mu
                  (b01 + b11*T + b21*pow(T, 2))*pow(abs(S), 0.5) +
                  (b02)*S;                                                         
    }

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

    //Apply boundary conditions
    w.ProjectCoefficient(w_coeff);
    w.ProjectBdrCoefficient(closed_down, ess_bdr_w);
 
    psi.ProjectCoefficient(psi_coeff);
    psi.ProjectBdrCoefficient(psi_in, bdr_psi_in);
    psi.ProjectBdrCoefficient(psi_out, bdr_psi_out);
    psi.ProjectBdrCoefficient(closed_down, bdr_psi_closed_down);
    psi.ProjectBdrCoefficient(closed_up, bdr_psi_closed_up);

    //Define the non-constant RHS
    if (f) delete f;
    f = new ParLinearForm(&fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(k_r_Theta_dr));
    f->AddDomainIntegrator(new DomainLFIntegrator(k_r_Phi_dr));
    f->Assemble();

    //Define non-constant bilinear forms of the system
    if (d) delete d;
    d = new ParBilinearForm (&fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(neg_Eta));
    d->AddDomainIntegrator(new ConvectionIntegrator(neg_Eta_r_inv_hat));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, psi, *f, Operator::DIAG_KEEP);
    d->Finalize();
    if (D) delete D;
    D = d->ParallelAssemble();

    if (ct) delete ct;
    ct = new ParMixedBilinearForm(&fespace, &fespace);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator);
    ct->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    ct->Assemble();
    ct->EliminateTrialDofs(ess_bdr_w, w, *f);
    ct->EliminateTestDofs(ess_bdr_psi);
    ct->Finalize();
    if (Ct) delete Ct;
    Ct = ct->ParallelAssemble();

    //Transfer to TrueDofs
    f->ParallelAssemble(B.GetBlock(1));
}
