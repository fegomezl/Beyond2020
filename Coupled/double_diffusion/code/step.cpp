#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(X);
    ode_solver->Step(X, t, dt);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        theta->Distribute(&(X.GetBlock(0)));
        phi->Distribute(&(X.GetBlock(1)));

        //Calculate phases
        double T_f;
        for (int ii = 0; ii < phase->Size(); ii++){
            T_f = T_fun((*phi)(ii));
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

void Conduction_Operator::SetParameters(const BlockVector &X){
    //Recover actual information
    theta.Distribute(&(X.GetBlock(0)));
    phi.Distribute(&(X.GetBlock(1)));

    //Associate the values of each auxiliar function
    double DT = 0.;
    for (int ii = 0; ii < phase.Size(); ii++){
        DT = theta(ii) - T_fun(phi(ii)); 
        theta(ii) = DT;
        phase(ii) = 0.5*(1 + tanh(5*config.invDeltaT*DT));
    }

    aux_A.Set(config.a_l-config.a_s, phase); aux_A += config.a_s;
    aux_D.Set(config.d_l-config.d_s, phase); aux_D += config.d_s;
    aux_L.Set(config.L_l-config.L_s, phase); aux_L += config.L_s;

    //Set the associated coefficients
    GridFunctionCoefficient coeff_A(&aux_A);
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

    //Construct final coefficients
    coeff_rA.SetBCoef(coeff_A);
    coeff_rD.SetBCoef(coeff_D); 
    coeff_rL.SetBCoef(LDeltaT); 

    //Create corresponding bilinear forms
    delete m_theta;
    delete M_theta;
    delete M_e_theta;
    delete M_0_theta;
    m_theta = new ParBilinearForm(&fespace);
    m_theta->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m_theta->AddDomainIntegrator(new MassIntegrator(coeff_rL));
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
    k_theta->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA));
    k_theta->Assemble();
    k_theta->Finalize();
    K_0_theta = k_theta->ParallelAssemble();

    delete k_phi;
    delete K_0_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->Assemble();
    k_phi->Finalize();
    K_0_phi = k_phi->ParallelAssemble();
}
