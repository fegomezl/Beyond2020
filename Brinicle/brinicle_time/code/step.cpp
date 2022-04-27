#include "header.h"

//Evolve the simulation one time step 
void Artic_sea::time_step(){

    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *rVelocity);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X, t);
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
            (*phase)(ii) = Phase((*temperature)(ii), (*salinity)(ii));
    
        //Adimentionalize stream function
        if (config.rescale){
            double stream_local_max = stream->Max(), stream_max;
            double stream_local_min = stream->Min(), stream_min;
            MPI_Allreduce(&stream_local_max, &stream_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&stream_local_min, &stream_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            for (int ii = 0; ii < stream->Size(); ii++)
                    (*stream)(ii) = ((*stream)(ii)-stream_min)/(stream_max-stream_min);
        }

        //Print fields
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

//Update of the solver on each iteration
void Transport_Operator::SetParameters(const BlockVector &X, const Vector &rVelocity){
    //Recover current information
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

//Update of the solver on each iteration
void Flow_Operator::SetParameters(const BlockVector &X, double t){

    //Recover current information
    temperature.SetFromTrueDofs(X.GetBlock(0));
    salinity.SetFromTrueDofs(X.GetBlock(1));

    //Calculate impermeability and density coefficients
    for (int ii = 0; ii < impermeability.Size(); ii++){
        double T = temperature(ii);
        double S = salinity(ii);

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
    ConstantCoefficient coeff_buoyancy_constant(constants.BuoyancyCoefficient);
    ProductCoefficient coeff_buoyancy(coeff_buoyancy_constant, coeff_density_dr);
    ProductCoefficient coeff_r_buoyancy(coeff_r, coeff_buoyancy);
    
    //Apply boundary conditions
    coeff_vorticity.SetTime(t);
    vorticity_boundary.ProjectCoefficient(coeff_vorticity);
    vorticity_boundary.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);

    coeff_stream.SetTime(t);
    coeff_stream_in.SetTime(t);
    coeff_stream_out.SetTime(t);
    coeff_stream_closed_down.SetTime(t);
    coeff_stream_closed_up.SetTime(t);
    stream_boundary.ProjectCoefficient(coeff_stream);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_in);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_out, ess_bdr_out);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_closed_down);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_closed_up);

    //Define non-constant RHS
    if (B0) delete B0;
    ParLinearForm b0(&fespace_H1);
    ConstantCoefficient Zero(0.);
    b0.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    b0.Assemble();

    if (B1) delete B1;
    ParLinearForm b1(&fespace_H1);
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy));
    b1.Assemble();

    //Define non-constant bilinear forms of the system
    if (A00) delete A00;
    ParBilinearForm a00(&fespace_H1);
    a00.AddDomainIntegrator(new MassIntegrator);
    a00.Assemble();
    a00.EliminateEssentialBC(ess_bdr_0, vorticity_boundary, b0, Operator::DIAG_ONE);
    a00.Finalize();
    A00 = a00.ParallelAssemble();

    if (A10) delete A10;
    ParMixedBilinearForm a10(&fespace_H1, &fespace_H1);
    a10.AddDomainIntegrator(new MixedGradGradIntegrator);
    a10.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a10.Assemble();
    a10.EliminateTrialDofs(ess_bdr_0, vorticity_boundary, b1);
    a10.EliminateTestDofs(ess_bdr_1);    
    a10.Finalize();
    A10 = a10.ParallelAssemble();

    if (A01) delete A01;
    ParMixedBilinearForm a01(&fespace_H1, &fespace_H1);
    a01.AddDomainIntegrator(new MixedGradGradIntegrator);
    a01.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a01.Assemble();
    a01.EliminateTrialDofs(ess_bdr_1, stream_boundary, b0);
    a01.EliminateTestDofs(ess_bdr_0);    
    a01.Finalize();
    A01 = a01.ParallelAssemble();

    if (A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_impermeability));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_impermeability_r_inv_hat));
    a11.Assemble();
    a11.EliminateEssentialBC(ess_bdr_1, stream_boundary, b1, Operator::DIAG_KEEP);
    a11.Finalize();
    A11 = a11.ParallelAssemble();

    //Transfer to TrueDofs
    B0 = b0.ParallelAssemble();
    B1 = b1.ParallelAssemble();
}
