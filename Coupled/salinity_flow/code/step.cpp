#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(*Temperature, *Salinity, *rVelocity);
    ode_solver->Step(*Salinity, t, dt);

    flow_oper->SetParameters(*Temperature, *Salinity);
    flow_oper->Solve(Y, *Velocity, *rVelocity);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_print++;

        //Update information
        temperature->Distribute(Temperature);
        salinity->Distribute(Salinity);
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

void Transport_Operator::SetParameters(const Vector &Temperature, const Vector &Salinity, const Vector &rVelocity){
    //Recover actual information
    temperature.SetFromTrueDofs(Temperature);
    salinity.SetFromTrueDofs(Salinity);
    rvelocity.SetFromTrueDofs(rVelocity); 

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < phase.Size(); ii++){
        salt_diffusivity(ii) = SaltDiffusivity(temperature(ii), salinity(ii));;

        phase(ii) = Phase(temperature(ii), salinity(ii));
        temperature(ii) = temperature(ii) - FusionPoint(salinity(ii));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_D(&salt_diffusivity);

    //Construct final coefficients
    coeff_rD.SetBCoef(coeff_D); 

    coeff_rV.SetGridFunction(&rvelocity);

    //Create transport matrix
    if (K) delete K;
    ParBilinearForm k(&fespace_H1);
    k.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k.AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k.Assemble();
    k.Finalize();
    K = k.ParallelAssemble();    
}

void Flow_Operator::SetParameters(const Vector &Temperature, const Vector &Salinity){
    //Update information
    temperature.SetFromTrueDofs(Temperature);
    salinity.SetFromTrueDofs(Salinity);

    temperature.GetDerivative(1, 0, temperature_dr);
    salinity.GetDerivative(1, 0, salinity_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < phase.Size(); ii++){
        double T = temperature(ii);
        double S = salinity(ii);

        impermeability(ii) = Impermeability(T, S);
        temperature(ii) = ExpansivityTemperature(T, S);
        salinity(ii) = ExpansivitySalinity(T, S);
    }

    //Properties coefficients
    GridFunctionCoefficient coeff_impermeability(&impermeability);
    ProductCoefficient coeff_neg_impermeability(-1., coeff_impermeability);

    GridFunctionCoefficient coeff_temperature_dr(&temperature_dr);
    GridFunctionCoefficient coeff_expansivity_temperature(&temperature);
    ProductCoefficient coeff_buoyancy_temperature(coeff_expansivity_temperature, coeff_temperature_dr);

    GridFunctionCoefficient coeff_salinity_dr(&salinity_dr);
    GridFunctionCoefficient coeff_expansivity_salinity(&salinity);
    ProductCoefficient coeff_buoyancy_salinity(coeff_expansivity_salinity, coeff_salinity_dr);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient coeff_neg_impermeability_r_inv_hat(coeff_neg_impermeability, coeff_r_inv_hat);
    ProductCoefficient coeff_r_buoyancy_temperature(coeff_r, coeff_buoyancy_temperature);
    ProductCoefficient coeff_r_buoyancy_salinity(coeff_r, coeff_buoyancy_salinity);
    
    //Apply boundary conditions
    vorticity_boundary.ProjectCoefficient(coeff_vorticity);
    vorticity_boundary.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);

    stream_boundary.ProjectCoefficient(coeff_stream);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_in);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_out, ess_bdr_out);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_closed_down);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_closed_up);

    //Define the non-constant RHS
    if (B1) delete B1;
    ParLinearForm b1(&fespace_H1);
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_temperature));
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_salinity));
    b1.Assemble();

    //Define non-constant bilinear forms of the system
    if (A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_impermeability));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_impermeability_r_inv_hat));
    a11.Assemble();
    a11.EliminateEssentialBC(ess_bdr_1, stream_boundary, b1, Operator::DIAG_KEEP);
    a11.Finalize();
    A11 = a11.ParallelAssemble();

    if (A10) delete A10;
    ParMixedBilinearForm a10(&fespace_H1, &fespace_H1);
    a10.AddDomainIntegrator(new MixedGradGradIntegrator);
    a10.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a10.Assemble();
    a10.EliminateTrialDofs(ess_bdr_0, vorticity_boundary, b1);
    a10.EliminateTestDofs(ess_bdr_1);    
    a10.Finalize();
    A10 = a10.ParallelAssemble();

    //Transfer to TrueDofs
    B1 = b1.ParallelAssemble();
}
