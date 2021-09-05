#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(*Theta, *V);
    ode_solver->Step(*Theta, t, dt);

    flow_oper->SetParameters(*Theta);
    flow_oper->Solve(*W, *Psi, *V);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        theta->Distribute(Theta);
        w->Distribute(W);
        psi->Distribute(Psi);
        v->Distribute(V);

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

void Conduction_Operator::SetParameters(const Vector &X, const Vector &V){
    //Read Velocity
    v.SetFromTrueDofs(V);
    rV.SetGridFunction(&v);

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
    coeff_rCLV.SetACoef(coeff_CL);  coeff_rCLV.SetBCoef(rV);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

    //Create corresponding bilinear forms
    if (m) delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    m->Assemble();
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    if (k) delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k->Assemble();
    k->Finalize();
}

void Flow_Operator::SetParameters(const Vector &Theta){
    //Update temperature coefficients
    theta.SetFromTrueDofs(Theta);
    theta.GetDerivative(1, 0, theta_dr);

    theta_eta.SetFromTrueDofs(Theta);
    for (int ii = 0; ii < theta_eta.Size(); ii++){
        theta_eta(ii) = 0.5*(1 + tanh(5*config.invDeltaT*(theta_eta(ii) - config.T_f)));
        theta_eta(ii) = config.epsilon_eta + pow(1-theta_eta(ii), 2)/(pow(theta_eta(ii), 3) + config.epsilon_eta);
    }

    //Properties coefficients
    GridFunctionCoefficient Theta_dr(&theta_dr);
    ProductCoefficient k_Theta_dr(30., Theta_dr);

    GridFunctionCoefficient eta(&theta_eta);
    ProductCoefficient neg_eta(-1., eta);

    //Rotational coupled coefficients
    ProductCoefficient k_r_Theta_dr(r, k_Theta_dr);
    ScalarVectorProductCoefficient neg_eta_r_inv_hat(neg_eta, r_inv_hat);

    //Apply boundary conditions
    w.ProjectCoefficient(w_coeff);
    psi.ProjectCoefficient(psi_coeff);

    //Define the RHS
    if (g) delete g;
    g = new ParLinearForm(&fespace);
    g->Assemble();

    if (f) delete f;
    f = new ParLinearForm(&fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(k_r_Theta_dr));
    f->Assemble();

    //Define bilinear forms of the system
    if (m) delete m;
    m = new ParBilinearForm (&fespace);
    m->AddDomainIntegrator(new MassIntegrator);
    m->Assemble();
    m->EliminateEssentialBC(ess_bdr_w, w, *g, Operator::DIAG_ONE);
    m->Finalize();
    M = m->ParallelAssemble();

    if (d) delete d;
    d = new ParBilinearForm (&fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(neg_eta));
    d->AddDomainIntegrator(new ConvectionIntegrator(neg_eta_r_inv_hat));
    d->Assemble();
    d->EliminateEssentialBC(ess_bdr_psi, psi, *f, Operator::DIAG_KEEP);
    d->Finalize();
    D = d->ParallelAssemble();

    if (c) delete c;
    c = new ParMixedBilinearForm (&fespace, &fespace);
    c->AddDomainIntegrator(new MixedGradGradIntegrator);
    c->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    c->Assemble();
    c->EliminateTrialDofs(ess_bdr_psi, psi, *g);
    c->EliminateTestDofs(ess_bdr_w);
    c->Finalize();
    C = c->ParallelAssemble();

    if (ct) delete ct;
    ct = new ParMixedBilinearForm(&fespace, &fespace);
    ct->AddDomainIntegrator(new MixedGradGradIntegrator);
    ct->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(r_inv_hat));
    ct->Assemble();
    ct->EliminateTrialDofs(ess_bdr_w, w, *f);
    ct->EliminateTestDofs(ess_bdr_psi);
    ct->Finalize();
    Ct = ct->ParallelAssemble();

    //Transfer to TrueDofs
    w.ParallelAssemble(X.GetBlock(0));
    psi.ParallelAssemble(X.GetBlock(1));

    g->ParallelAssemble(B.GetBlock(0));
    f->ParallelAssemble(B.GetBlock(1));
}
