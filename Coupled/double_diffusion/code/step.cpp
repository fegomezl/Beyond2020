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

void Conduction_Operator::SetParameters(const BlockVector &X){
    //Recover actual information
    aux_theta.Distribute(&(X.GetBlock(0)));
    aux_phi.Distribute(&(X.GetBlock(1)));

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux_theta.Size(); ii++){
        aux_theta(ii) -= config.T_f + T_fun(aux_phi(ii));

        if (aux_theta(ii) > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
            aux_D(ii) = config.D_l;
        } else {
            aux_C(ii) = config.c_s + delta_c_s_fun(aux_theta(ii), aux_phi(ii));
            aux_K(ii) = config.k_s;
            aux_D(ii) = config.D_s;
        }

        aux_L(ii) = 0.5*config.L*(1 + tanh(5*config.invDeltaT*aux_theta(ii)));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);

    //Construct latent heat term
    GradientGridFunctionCoefficient dT(&aux_theta);
    GradientGridFunctionCoefficient dH(&aux_L);
    
    dHdT.SetACoef(dH);  dT_2.SetACoef(dT);
    dHdT.SetBCoef(dT);  dT_2.SetBCoef(dT);

    SumCoefficient dT_2e(config.EpsilonT, dT_2);
    
    PowerCoefficient inv_dT_2(dT_2e, -1.);
    ProductCoefficient LDeltaT(dHdT, inv_dT_2);

    SumCoefficient coeff_CL(coeff_C, LDeltaT);

    //Construct final coefficients
    coeff_rC.SetBCoef(coeff_CL);
    coeff_rK.SetBCoef(coeff_K); 
    coeff_rD.SetBCoef(coeff_D); 

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
    k_theta->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_theta), robin_bdr_theta);
    k_theta->Assemble();
    k_theta->Finalize();
    K_0_theta = k_theta->ParallelAssemble();

    if(k_phi) delete k_phi;
    if(K_0_phi) delete K_0_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->AddBoundaryIntegrator(new MassIntegrator(r_robin_h_phi), robin_bdr_phi);
    k_phi->Assemble();
    k_phi->Finalize();
    K_0_phi = k_phi->ParallelAssemble();
}
