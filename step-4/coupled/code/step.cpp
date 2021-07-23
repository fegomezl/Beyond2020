#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    oper_T->SetParameters(X);
    ode_solver->Step(X, t, dt);

    //Solve the flow problem
    theta->Distribute(&(X.GetBlock(0)));
    Theta = theta->GetTrueDofs();

    flow_oper->Solve(Theta);
    (*x_psi) = flow_oper->GetStream();
    (*x_v) = flow_oper->GetVelocity();
    (*x_w) = flow_oper->GetVorticity();

    oper_T->UpdateVelocity(x_v);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = config.T_f + T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
        }

        //Graph
        phi->Distribute(&(X.GetBlock(1)));
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
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);

    GridFunctionCoefficient coeff_L(&aux_theta);
    SumCoefficient coeff_CL(coeff_C, coeff_L);
    coeff_rCL.SetBCoef(coeff_CL);

    coeff_rK.SetBCoef(coeff_K); dt_coeff_rK.SetBCoef(coeff_rK);

    coeff_rD.SetBCoef(coeff_D); dt_coeff_rD.SetBCoef(coeff_rD);

    coeff_rCLV.SetACoef(coeff_CL);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

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
    k_theta->AddBoundaryIntegrator(new MassIntegrator(neg_r_robin_h_theta), robin_bdr_theta);
    k_theta->Assemble();
    k_theta->Finalize();

    delete k_phi;
    k_phi = new ParBilinearForm(&fespace);
    k_phi->AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k_phi->AddBoundaryIntegrator(new MassIntegrator(neg_r_robin_h_phi), robin_bdr_phi);
    k_phi->Assemble();
    k_phi->Finalize();
}

void Conduction_Operator::UpdateVelocity(const ParGridFunction *v_flow){
    v = *v_flow;
    coeff_rV.SetGridFunction(&v);
    dt_coeff_rV.SetBCoef(coeff_rV);
    coeff_rCLV.SetBCoef(coeff_rV);
}
