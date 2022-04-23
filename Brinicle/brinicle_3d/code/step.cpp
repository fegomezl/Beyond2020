#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *rVelocity);
    ode_solver->Step(X, t, dt);

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
    m0.AddDomainIntegrator(new MassIntegrator(coeff_M));
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
    k0.AddDomainIntegrator(new DiffusionIntegrator(coeff_D0));
    //k0.AddDomainIntegrator(new ConvectionIntegrator(coeff_rMV));
    k0.Assemble();
    k0.Finalize();
    K0 = k0.ParallelAssemble();    

    if (K1) delete K1;
    ParBilinearForm k1(&fespace_H1);
    k1.AddDomainIntegrator(new DiffusionIntegrator(coeff_D1));
    //k1.AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k1.Assemble();
    k1.Finalize();
    K1 = k1.ParallelAssemble();
}
