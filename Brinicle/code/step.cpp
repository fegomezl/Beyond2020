#include "header.h"

//Evolve the simulation one time step 
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

    //Print visualization on certain steps
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
            (*phase)(ii) = FusionPoint((*temperature)(ii), (*salinity)(ii));

        //Print fields
        paraview_out->SetCycle(vis_print);
        paraview_out->SetTime(t*t_ref);
        paraview_out->Save();
    }

    //Print the system state
    double percentage = 100*t/config.t_final;
    string progress = to_string((int)percentage)+"%";
    if (config.master){
        cout.flush();
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << vis_print << setw(12)
             << dt*t_ref << setw(12)
             << t*t_ref  << setw(12)
             << progress << "\r";
    }
}

//Update of the solver on each iteration
void Transport_Operator::SetParameters(const BlockVector &X, const Vector &rVelocity){
    //Recover current information
    temperature.SetFromTrueDofs(X.GetBlock(0));
    salinity.SetFromTrueDofs(X.GetBlock(1));
    rvelocity.SetFromTrueDofs(rVelocity); 

    //Associate the values of each auxiliar function
    double T, S = 0.;
    for (int ii = 0; ii < phase.Size(); ii++){
        T = temperature(ii);
        S = salinity(ii);

        heat_diffusivity(ii) = HeatDiffusivity(T, S);
        salt_diffusivity(ii) = SaltDiffusivity(T, S);
        phase(ii)            = Phase(T, S);
        temperature(ii)      = FusionPoint(T, S);
    }

    //Set the associated coefficients
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

    ProductCoefficient coeff_latent(constants.Stefan, DeltaT);
    SumCoefficient coeff_M(1., coeff_latent);

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

//Update of the solver on each iteration
void Flow_Operator::SetParameters(const BlockVector &X){

    //Recover current information
    temperature.SetFromTrueDofs(X.GetBlock(0));
    salinity.SetFromTrueDofs(X.GetBlock(1));

    //Calculate impermeability and density coefficients
    double T, S = 0.;
    for (int ii = 0; ii < density.Size(); ii++){
        T = temperature(ii);
        S = salinity(ii);

        impermeability(ii) = Impermeability(T, S);
        density(ii) = Density(T, S);
    }
    
    //Calculate gradient of the density field
    density.GetDerivative(1, 0, density_dr);

    //Create impermeability and buoyancy coefficients
    GridFunctionCoefficient coeff_impermeability(&impermeability);
    ProductCoefficient coeff_neg_impermeability(-1., coeff_impermeability);
    ScalarVectorProductCoefficient coeff_neg_impermeability_r_inv_hat(coeff_neg_impermeability, coeff_r_inv_hat);

    GridFunctionCoefficient coeff_density_dr(&density_dr);
    ProductCoefficient coeff_buoyancy(constants.BuoyancyCoefficient, coeff_density_dr);
    ProductCoefficient coeff_r_buoyancy(coeff_r, coeff_buoyancy);
    
    //Apply boundary conditions
    vorticity.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);

    stream.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_in);
    stream.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_closed_down);
    stream.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_closed_up);

    //Define non-constant RHS
    if (B1) delete B1;
    ParLinearForm b1(&fespace_H1);
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy));
    b1.Assemble();

    //Define non-constant bilinear forms of the system
    if (A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_impermeability));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_impermeability_r_inv_hat));
    a11.Assemble();
    a11.EliminateEssentialBC(ess_bdr_1, stream, b1, Operator::DIAG_KEEP);
    a11.Finalize();
    A11 = a11.ParallelAssemble();

    if (A10) delete A10;
    ParMixedBilinearForm a10(&fespace_H1, &fespace_H1);
    a10.AddDomainIntegrator(new MixedGradGradIntegrator);
    a10.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a10.Assemble();
    a10.EliminateTrialDofs(ess_bdr_0, vorticity, b1);
    a10.EliminateTestDofs(ess_bdr_1);    
    a10.Finalize();
    A10 = a10.ParallelAssemble();

    //Transfer to TrueDofs
    B1 = b1.ParallelAssemble();
}
