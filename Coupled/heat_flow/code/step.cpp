#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(*Theta, *rV);
    ode_solver->Step(*Theta, t, dt);

    flow_oper->SetParameters(*Theta);
    flow_oper->Solve(*W, *Psi, *V, *rV);

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

void Conduction_Operator::SetParameters(const Vector &X, const Vector &rV){
    //Read Velocity
    rv.SetFromTrueDofs(rV);
    coeff_rV.SetGridFunction(&rv);

    //Create the auxiliar grid functions
    aux.SetFromTrueDofs(X);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        aux(ii) -= config.T_f;

        if (aux(ii) > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
        } else {
            aux_C(ii) = config.c_s;
            aux_K(ii) = config.k_s;
        }

        aux_L(ii) = 0.5*config.L*(1 + tanh(5*config.invDeltaT*aux(ii)));
    }

    //Set the grid coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);

    //Construct latent heat term
    GradientGridFunctionCoefficient dT(&aux);
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
    coeff_rCV.SetACoef(coeff_CL);   coeff_rCV.SetBCoef(coeff_rV);   

    //Create corresponding bilinear forms
    if (m) delete m;
    if (M) delete M;
    if (M_e) delete M_e;
    if (M_0) delete M_0;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();
    M_e = M->EliminateRowsCols(ess_tdof_list);
    M_0 = m->ParallelAssemble();

    M_prec.SetOperator(*M);
    M_solver.SetOperator(*M);

    if (k) delete k;
    if (K_0) delete K_0;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCV));
    k->Assemble();
    k->Finalize();
    K_0 = k->ParallelAssemble();    
}

void Flow_Operator::SetParameters(const Vector &Theta){
    //Update temperature coefficients
    theta.SetFromTrueDofs(Theta);
    theta.GetDerivative(1, 0, theta_dr);
    for (int ii = 0; ii < theta.Size(); ii++)
        theta(ii) = 3.8*theta(ii) - 16;     // gb/u

    theta_eta.SetFromTrueDofs(Theta);
    for (int ii = 0; ii < theta_eta.Size(); ii++){
        theta_eta(ii) = 0.5*(1 + tanh(5*config.invDeltaT*(theta_eta(ii) - config.T_f)));
        theta_eta(ii) = config.EpsilonEta + pow(1-theta_eta(ii), 2)/(pow(theta_eta(ii), 3) + config.EpsilonEta);
    }

    //Properties coefficients
    GridFunctionCoefficient Theta_dr(&theta_dr);
    GridFunctionCoefficient k(&theta);
    ProductCoefficient k_Theta_dr(k, Theta_dr);

    GridFunctionCoefficient eta(&theta_eta);
    ProductCoefficient neg_eta(-1., eta);

    //Rotational coupled coefficients
    ProductCoefficient k_r_Theta_dr(coeff_r, k_Theta_dr);
    ScalarVectorProductCoefficient neg_eta_r_inv_hat(neg_eta, r_inv_hat);

    //Apply boundary conditions
    w.ProjectCoefficient(w_coeff);
    psi.ProjectCoefficient(psi_coeff);

    //Define the non-constant RHS
    if (f) delete f;
    f = new ParLinearForm(&fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(k_r_Theta_dr));
    f->Assemble();

    //Define non-constant bilinear forms of the system
    if (d) delete d;
    d = new ParBilinearForm (&fespace);
    d->AddDomainIntegrator(new DiffusionIntegrator(neg_eta));
    d->AddDomainIntegrator(new ConvectionIntegrator(neg_eta_r_inv_hat));
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
