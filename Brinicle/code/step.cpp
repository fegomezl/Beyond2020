#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(X, *rV);
    ode_solver->Step(X, t, dt);

    /*//Normalize the salinity
    double m_in, m_out;
    ConstantCoefficient low_cap(0.);
    ParGridFunction phi_aux(*phi);
    phi_aux.SetFromTrueDofs(X.GetBlock(1));
    phi_aux -= phi_in;

    m_in = phi_aux.ComputeL1Error(low_cap, irs);
    for (int ii = 0; ii < phi_aux.Size(); ii++){
        if (phi_aux(ii) < 0)
            phi_aux(ii) = 0;
    }
    m_out = phi_aux.ComputeL1Error(low_cap, irs);
    phi_aux *= 2-m_in/m_out; 
    
    m_in = phi_aux.ComputeL1Error(low_cap, irs);
    for (int ii = 0; ii < phi_aux.Size(); ii++){
        if (phi_aux(ii) > phi_out-phi_in)
            phi_aux(ii) = phi_out-phi_in;
    }
    m_out = phi_aux.ComputeL1Error(low_cap, irs);
    phi_aux *= m_in/m_out; 

    phi_aux += phi_in;
    phi_aux.GetTrueDofs(X.GetBlock(1));*/

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
            double T_f = config.T_f + T_fun((*phi)(ii));
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
    theta.Distribute(X.GetBlock(0));
    phi.Distribute(X.GetBlock(1));
    rv.Distribute(rV); 

    //Associate the values of each auxiliar function
    double DT = 0.;
    for (int ii = 0; ii < phase.Size(); ii++){
        DT = theta(ii) - config.T_f - T_fun(phi(ii));
        theta(ii) = DT;
        phase(ii) = 0.5*(1 + tanh(5*config.invDeltaT*DT));
    }

    aux_C.Set(config.c_l-config.c_s, phase); aux_C += config.c_s;
    aux_K.Set(config.k_l-config.k_s, phase); aux_K += config.k_s;
    aux_D.Set(config.d_l-config.d_s, phase); aux_D += config.d_s;
    aux_L.Set(config.L_l-config.L_s, phase); aux_L += config.L_s;

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
    delete m_theta;
    delete M_theta;
    delete M_e_theta;
    delete M_0_theta;
    m_theta = new ParBilinearForm(&fespace);
    m_theta->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m_theta->Assemble();
    m_theta->Finalize();
    M_theta = m_theta->ParallelAssemble();
    M_e_theta = M_theta->EliminateRowsCols(ess_tdof_theta);
    M_0_theta = m_theta->ParallelAssemble();

    M_theta_prec.SetOperator(*M_theta);
    M_theta_solver.SetOperator(*M_theta);

    delete m_phi;
    delete M_phi;
    delete M_e_phi;
    delete M_0_phi;
    m_phi = new ParBilinearForm(&fespace);
    m_phi->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m_phi->Assemble();
    m_phi->Finalize();
    M_phi = m_phi->ParallelAssemble();
    M_e_phi = M_phi->EliminateRowsCols(ess_tdof_phi);
    M_0_phi = m_phi->ParallelAssemble();

    M_phi_prec.SetOperator(*M_phi);
    M_phi_solver.SetOperator(*M_phi);

    delete k_theta;
    delete K_0_theta;
    k_theta = new ParBilinearForm(&fespace);
    k_theta->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k_theta->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCV));
    k_theta->Assemble();
    k_theta->Finalize();
    K_0_theta = k_theta->ParallelAssemble();

    delete k_phi;
    delete K_0_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k_phi->Assemble();
    k_phi->Finalize();
    K_0_phi = k_phi->ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Update information
    theta.Distribute(X.GetBlock(0));
    phi.Distribute(X.GetBlock(1));

    theta.GetDerivative(1, 0, theta_dr);
    phi.GetDerivative(1, 0, phi_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < phase.Size(); ii++){
        double T = theta(ii);
        double S = phi(ii);
        double P = 0.5*(1 + tanh(5*config.invDeltaT*(theta(ii) - config.T_f - T_fun(phi(ii))))); 

        eta(ii) = config.EpsilonEta + pow(1-P, 2)/(pow(P, 3) + config.EpsilonEta);

        theta(ii) = delta_rho_t_fun(T, S);
        phi(ii) = delta_rho_p_fun(T, S);
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
    w.ProjectBdrCoefficient(closed_down, ess_bdr_w);
    w.ParallelProject(W);
 
    psi.ProjectBdrCoefficient(psi_in, bdr_psi_in);
    psi.ProjectBdrCoefficient(psi_out, bdr_psi_out);
    psi.ProjectBdrCoefficient(closed_down, bdr_psi_closed_down);
    psi.ProjectBdrCoefficient(closed_up, bdr_psi_closed_up);
    psi.ParallelProject(Psi);

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
