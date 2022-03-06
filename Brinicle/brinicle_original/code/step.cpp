#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *rVelocity);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *rVelocity);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_print++;

        //Update information
        temperature->Distribute(X.GetBlock(0));
        salinity->Distribute(X.GetBlock(1));
        vorticity->Distribute(Y.GetBlock(0));
        stream->Distribute(Y.GetBlock(1));
        velocity->Distribute(Velocity);
        rvelocity->Distribute(rVelocity);
        
        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++)
            (*phase)(ii) = Phase((*temperature)(ii), (*salinity)(ii));
    
        //Normalize stream
        if (config.rescale){
            double stream_local_max = stream->Max(), stream_max;
            double stream_local_min = stream->Min(), stream_min;
            MPI_Allreduce(&stream_local_max, &stream_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&stream_local_min, &stream_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            for (int ii = 0; ii < stream->Size(); ii++)
                    (*stream)(ii) = ((*stream)(ii)-stream_min)/(stream_max-stream_min);
        }

        //Graph
        paraview_out->SetCycle(vis_print);
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

void Transport_Operator::SetParameters(const BlockVector &X, const Vector &rVelocity){
    //Recover actual information
    temperature.SetFromTrueDofs(X.GetBlock(0));
    salinity.SetFromTrueDofs(X.GetBlock(1));
    rvelocity.SetFromTrueDofs(rVelocity); 

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < phase.Size(); ii++){
        heat_inertia(ii) = HeatInertia(temperature(ii), salinity(ii));
        heat_diffusivity(ii) = HeatDiffusivity(temperature(ii), salinity(ii));
        salt_diffusivity(ii) = SaltDiffusivity(temperature(ii), salinity(ii));;

        phase(ii) = Phase(temperature(ii), salinity(ii));
        temperature(ii) = temperature(ii) - FusionPoint(salinity(ii));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_I(&heat_inertia);
    GridFunctionCoefficient coeff_D0(&heat_diffusivity);
    GridFunctionCoefficient coeff_D1(&salt_diffusivity);

    //Construct latent heat term
    GradientGridFunctionCoefficient coeff_dT(&temperature);
    GradientGridFunctionCoefficient coeff_dP(&phase);
    
    coeff_dPdT.SetACoef(coeff_dP);  coeff_dT_2.SetACoef(coeff_dT);
    coeff_dPdT.SetBCoef(coeff_dT);  coeff_dT_2.SetBCoef(coeff_dT);

    SumCoefficient coeff_dT_2e(Epsilon, coeff_dT_2);
    
    PowerCoefficient coeff_inv_dT_2(coeff_dT_2e, -1.);
    ProductCoefficient DeltaT(coeff_dPdT, coeff_inv_dT_2);

    SumCoefficient coeff_M(coeff_I, DeltaT);

    //Construct final coefficients
    coeff_rM.SetBCoef(coeff_M);
    coeff_rD0.SetBCoef(coeff_D0); 
    coeff_rD1.SetBCoef(coeff_D1); 

    coeff_rV.SetGridFunction(&rvelocity);
    coeff_rMV.SetACoef(coeff_M);
    coeff_rMV.SetBCoef(coeff_rV);

    //Create corresponding bilinear forms
    if (M0) delete M0;
    if (M0_e) delete M0_e;
    if (M0_o) delete M0_o;
    ParBilinearForm m0(&fespace_H1);
    m0.AddDomainIntegrator(new MassIntegrator(coeff_rM));
    m0.Assemble();
    m0.Finalize();
    M0 = m0.ParallelAssemble();
    M0_e = M0->EliminateRowsCols(ess_tdof_0);
    M0_o = m0.ParallelAssemble();

    M0_prec.SetOperator(*M0);
    M0_solver.SetOperator(*M0);

    //Create transport matrix
    if (K0) delete K0;
    ParBilinearForm k0(&fespace_H1);
    k0.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD0));
    k0.AddDomainIntegrator(new ConvectionIntegrator(coeff_rMV));
    k0.Assemble();
    k0.Finalize();
    K0 = k0.ParallelAssemble();    

    if (K1) delete K1;
    ParBilinearForm k1(&fespace_H1);
    k1.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD1));
    k1.AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k1.Assemble();
    k1.Finalize();
    K1 = k1.ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Update information
    theta.SetFromTrueDofs(X.GetBlock(0));
    phi.SetFromTrueDofs(X.GetBlock(1));

    theta.GetDerivative(1, 0, theta_dr);
    phi.GetDerivative(1, 0, phi_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < phase.Size(); ii++){
        double T = theta(ii);
        double S = phi(ii);
        double P = 0.5*(1 + tanh(5*EpsilonInv*(theta(ii) - FusionPoint(phi(ii))))); 

        eta(ii) = Epsilon + pow(1-P, 2)/(pow(P, 3) + Epsilon);

        theta(ii) = ExpansivityTemperature(T, S);
        phi(ii) = ExpansivitySalinity(T, S);
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
