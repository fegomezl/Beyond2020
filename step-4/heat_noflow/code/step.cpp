#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(*Theta);
    ode_solver->Step(*Theta, t, dt);

    //Psi = (flow_oper->psi)->GetTrueDofs();
    //cond_oper->UpdateVelocity(*Psi, flow_oper->v);

    //Solve the flow problem
    //flow_oper->Solve(Theta);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Graph
        theta->SetFromTrueDofs(*Theta);
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

void Conduction_Operator::SetParameters(const Vector &X){
    //Create the auxiliar grid functions
    Vector Aux(X);
    if (config.T_f != 0)
        Aux -= config.T_f;
    aux.SetFromTrueDofs(Aux);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
        } else {
            aux_C(ii) = config.c_s;
            aux_K(ii) = config.k_s;
        }

        aux(ii) = config.L*config.invDeltaT*exp(-M_PI*pow(config.invDeltaT*aux(ii), 2));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_L(&aux);

    SumCoefficient coeff_CL(coeff_C, coeff_L);

    coeff_rK.SetBCoef(coeff_K); dt_coeff_rK.SetBCoef(coeff_rK); 
    coeff_rCL.SetBCoef(coeff_CL);
    coeff_rCLV.SetACoef(coeff_CL);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

    //Create corresponding bilinear forms
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    m->Assemble();
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k->Assemble();
    k->Finalize();
}

void Conduction_Operator::UpdateVelocity(const HypreParVector &Psi, ParGridFunction *v){
    psi.SetFromTrueDofs(Psi);
    gradpsi.SetGridFunction(&psi);
    rV.SetBCoef(gradpsi);
    v->Randomize();
    v->ProjectDiscCoefficient(rV, GridFunction::ARITHMETIC);
    coeff_rCLV.SetBCoef(rV);
}

void Flow_Operator::Update_T(const HypreParVector *Theta){
    //Update temperature coefficients
    if (theta_aux) delete theta_aux;
    theta_aux = new ParGridFunction(&fespace);
    theta_aux->SetFromTrueDofs(*Theta);
  
    if (theta_eta) delete theta_eta;
    theta_eta = new ParGridFunction(&fespace);
    theta_eta->SetFromTrueDofs(*Theta);
    for (int ii = 0; ii < theta_eta->Size(); ii++){
        (*theta_eta)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta_eta)(ii) - config.T_f)));
        (*theta_eta)(ii) = config.epsilon_eta + pow(1-(*theta_eta)(ii), 2)/(pow((*theta_eta)(ii), 3) + config.epsilon_eta);
    }
    
    if (theta_rho) delete theta_rho;
    theta_rho = new ParGridFunction(&fespace);
    theta_rho->SetFromTrueDofs(*Theta);
    for (int ii = 0; ii < theta_rho->Size(); ii++){
        if ((*theta_rho)(ii) > config.T_f)
            (*theta_rho)(ii) = 0.044*pow((*theta_rho)(ii), 2)
                             - 5.265*(*theta_rho)(ii)
                             + 19.886;
        else
            (*theta_rho)(ii) = 0.;
    }

    //Define local coefficients
    ConstantCoefficient neg(-1.);
  
    //Properties coefficients
    GridFunctionCoefficient eta(theta_eta);
    ScalarVectorProductCoefficient eta_r_inv_hat(eta, r_inv_hat);
  
    //Rotational coupled boundary coefficients
    ScalarVectorProductCoefficient neg_w_grad(neg, w_grad);
    ScalarVectorProductCoefficient neg_psi_grad(neg, psi_grad);
    InnerProductCoefficient r_inv_hat_w_grad(r_inv_hat, w_grad);

    InnerProductCoefficient eta_r_inv_hat_psi_grad(eta_r_inv_hat, psi_grad);
    ScalarVectorProductCoefficient eta_psi_grad(eta, psi_grad);
    ScalarVectorProductCoefficient neg_eta_psi_grad(eta, neg_psi_grad);
  
    //RHS coefficients
    GridFunctionCoefficient grad_rho(theta_rho);
    GradientGridFunctionCoefficient theta_grad(theta_aux);
    InnerProductCoefficient r_hat_theta_grad(r_hat, theta_grad);
    ProductCoefficient grad_rho_theta(grad_rho, r_hat_theta_grad);
    ProductCoefficient rF(r, grad_rho_theta);
    ProductCoefficient neg_rF(neg, rF);
  
    if(f) delete f;
    f = new ParLinearForm(&fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(neg_rF));
    f->AddDomainIntegrator(new DomainLFIntegrator(r_inv_hat_w_grad));
    f->AddDomainIntegrator(new DomainLFIntegrator(eta_r_inv_hat_psi_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(w_grad));
    f->AddDomainIntegrator(new DomainLFGradIntegrator(eta_psi_grad));
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_eta_psi_grad), ess_bdr_psi);
    f->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(neg_w_grad), ess_bdr_w);
    f->Assemble();
    f->ParallelAssemble(B.GetBlock(1));
  
    if(d) delete d;
    d = new ParBilinearForm(&fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(eta));
    d->AddDomainIntegrator(new ConvectionIntegrator(eta_r_inv_hat));
    d->Assemble();
    d->EliminateEssentialBCFromDofs(ess_tdof_list_psi);
    d->Finalize();
    D = d->ParallelAssemble();
}
