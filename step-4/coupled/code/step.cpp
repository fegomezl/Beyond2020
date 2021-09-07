#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(X, *V);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Z, *V);

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
        
        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = config.T_f + T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
        }

        //Graph
        paraview_out->SetCycle(vis_impressions);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }

    //Print the system state
    double percentage = 100*t/config.t_final;
    string progress = to_string((int)percentage)+"%";
    if (config.master){
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << dt << setw(12)
             << t  << setw(12)
             << progress << "\r";
        cout.flush();
    }
}

void Conduction_Operator::SetParameters(const BlockVector &X, const Vector &V){
    //Recover actual information
    aux_theta.SetFromTrueDofs(X.GetBlock(0));
    aux_phi.SetFromTrueDofs(X.GetBlock(1));
    v.SetFromTrueDofs(V); 

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux_phi.Size(); ii++){
        double T_f = config.T_f + T_fun(aux_phi(ii));
        if (aux_theta(ii) > T_f){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
            aux_D(ii) = config.D_l;
        } else {
            aux_C(ii) = config.c_s + delta_c_s_fun(aux_theta(ii), aux_phi(ii));
            aux_K(ii) = config.k_s;
            aux_D(ii) = config.D_s;
        }

        aux_theta(ii) = config.L*config.invDeltaT*exp(-M_PI*pow(config.invDeltaT*(aux_theta(ii) - T_f), 2));
    }


    //Set the associated coefficients
    coeff_rV.SetGridFunction(&v);
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);

    GridFunctionCoefficient coeff_L(&aux_theta);
    SumCoefficient coeff_CL(coeff_C, coeff_L);
    coeff_rCL.SetBCoef(coeff_CL);

    coeff_rK.SetBCoef(coeff_K); dt_coeff_rK.SetBCoef(coeff_rK);
    coeff_rCLV.SetACoef(coeff_CL); coeff_rCLV.SetBCoef(coeff_rV);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

    coeff_rD.SetBCoef(coeff_D); dt_coeff_rD.SetBCoef(coeff_rD);
    dt_coeff_rV.SetBCoef(coeff_rV);

    //Create corresponding bilinear forms
    delete m_theta;
    m_theta = new ParBilinearForm(&fespace);
    m_theta->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    m_theta->Assemble();
    m_theta->FormSystemMatrix(ess_tdof_list_theta, M_theta);
    M_theta_solver.SetOperator(M_theta);

    delete m_phi;
    m_phi = new ParBilinearForm(&fespace);
    m_phi->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m_phi->Assemble();
    m_phi->FormSystemMatrix(ess_tdof_list_phi, M_phi);
    M_phi_solver.SetOperator(M_phi);

    delete k_theta;
    k_theta = new ParBilinearForm(&fespace);
    k_theta->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k_theta->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k_theta->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_theta), robin_bdr_theta);
    k_theta->Assemble();
    k_theta->Finalize();

    delete k_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k_phi->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_phi), robin_bdr_phi);
    k_phi->Assemble();
    k_phi->Finalize();
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
        eta(ii) = config.epsilon_eta + pow(1-eta(ii), 2)/(pow(eta(ii), 3) + config.epsilon_eta);

        double a00 = 10.27542, a01 = -0.83195,
               a10 = -0.38667, a11 = 0.02801,
               a20 = 0.00624,  
               a30 = -0.00006;

        theta(ii) = (a00 + a10*T + a20*pow(T, 2) + a30*pow(T, 3))*S +               // k_t = g*b_t/mu
                    (a01 + a11*T)*pow(S, 1.5);                                  

        double b00 = -2061.94264, b01 = 68.54083, b02 = -24.43573,
               b10 = 10.27542,    b11 = -1.24793,
               b20 = -0.19333,    b21 = 0.02101,
               b30 = 0.00208,
               b40 = -0.00001;

        phi(ii) = (b00 + b10*T + b20*pow(T, 2) + b30*pow(T, 3) + b40*pow(T, 4)) +   // k_s = g*b_s/mu
                  (b01 + b11*T + b21*pow(T, 2))*pow(S, 0.5) +
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
    ProductCoefficient k_r_Theta_dr(r, k_Theta_dr);
    ProductCoefficient k_r_Phi_dr(r, k_Phi_dr);

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
