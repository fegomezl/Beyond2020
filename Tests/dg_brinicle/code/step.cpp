#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *r_Velocity);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        temperature->Distribute(X.GetBlock(0));
        salinity->Distribute(X.GetBlock(1));
        vorticity->Distribute(Y.GetBlock(0));
        stream->Distribute(Y.GetBlock(1));
        velocity->Distribute(Velocity);
        r_velocity->Distribute(r_Velocity);

        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++)
            (*phase)(ii) = Phase((*temperature)(ii), (*salinity)(ii));

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

void Transport_Operator::SetParameters(const BlockVector &X, const Vector &r_Velocity){
    //Recover actual information
    ParGridFunction temperature(&fespace_L2);
    ParGridFunction salinity(&fespace_L2);
    ParGridFunction r_velocity(&fespace_ND);

    temperature.Distribute(X.GetBlock(0));
    salinity.Distribute(X.GetBlock(1));
    r_velocity.Distribute(r_Velocity);

    //Create the auxiliar grid functions
    //ParGridFunction phase(&fespace_L2);
    //ParGridFunction heat_inertia(&fespace_L2);
    ParGridFunction heat_diffusivity(&fespace_L2);
    ParGridFunction salt_diffusivity(&fespace_L2);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < temperature.Size(); ii++){
        //heat_inertia(ii) = HeatInertia(temperature(ii), salinity(ii));
        heat_diffusivity(ii) = HeatDiffusivity(temperature(ii), salinity(ii));
        salt_diffusivity(ii) = SaltDiffusivity(temperature(ii), salinity(ii));
        //phase(ii) = Phase(temperature(ii), salinity(ii));
        //temperature(ii) -= FusionPoint(salinity(ii));
    }

    //temperature.ExchangeFaceNbrData();
    //phase.ExchangeFaceNbrData();
    //heat_inertia.ExchangeFaceNbrData();
    heat_diffusivity.ExchangeFaceNbrData();
    salt_diffusivity.ExchangeFaceNbrData();

    //Construct latent heat term
    //GradientGridFunctionCoefficient coeff_dT(&temperature), coeff_dH(&phase);
    //InnerProductCoefficient coeff_dHdT(coeff_dH, coeff_dT), coeff_dT_2(coeff_dT, coeff_dT);

    //SumCoefficient coeff_dT_2e(Epsilon, coeff_dT_2);
    
    //PowerCoefficient coeff_inv_dT_2(coeff_dT_2e, -1.);
    //ProductCoefficient coeff_delta(coeff_dHdT, coeff_inv_dT_2);

    //Set the grid coefficients
    //GridFunctionCoefficient coeff_M(&heat_inertia);
    //SumCoefficient coeff_MDelta(coeff_M, coeff_delta);

    GridFunctionCoefficient coeff_A0(&heat_diffusivity);
    GridFunctionCoefficient coeff_A1(&salt_diffusivity);

    ProductCoefficient coeff_rA0(coeff_r, coeff_A0);
    ProductCoefficient coeff_rA1(coeff_r, coeff_A1);
    //ProductCoefficient coeff_rMDelta(coeff_r, coeff_MDelta);

    VectorGridFunctionCoefficient coeff_rV(&r_velocity);
    //ScalarVectorProductCoefficient coeff_rMDeltaV(coeff_MDelta, coeff_rV);

    //Create mass matrix (Heat equation)
    m0 = new ParBilinearForm(&fespace_L2);
    m0->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m0->Assemble();
    m0->Finalize();
    M0 = m0->ParallelAssemble();
    M0_prec.SetOperator(*M0);
    M0_solver.SetOperator(*M0);

    //Create transport matrix
    config.sigma = 1.;
    config.kappa = 1.;
    config.eta = 0.;

    if (k0) delete k0;
    if (K0) delete K0;
    k0 = new ParBilinearForm(&fespace_L2);
    k0->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA0));
    k0->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k0->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rA0, config.sigma, config.kappa));
    k0->AddInteriorFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.));
    k0->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA0, config.sigma, config.kappa), ess_bdr_0);
    k0->AddBdrFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.), ess_bdr_0);
    k0->Assemble();
    k0->Finalize();
    K0 = k0->ParallelAssemble();    

    if (k1) delete k1;
    if (K1) delete K1;
    k1 = new ParBilinearForm(&fespace_L2);
    k1->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA1));
    k1->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k1->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rA1, config.sigma, config.kappa));
    k1->AddInteriorFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.));
    k1->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA1, config.sigma, config.kappa), ess_bdr_1);
    k1->AddBdrFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.), ess_bdr_1);
    k1->Assemble();
    k1->Finalize();
    K1 = k1->ParallelAssemble();    

    //Create RHS
    if (b0) delete b0;
    if (B0) delete B0;
    b0 = new ParLinearForm(&fespace_L2);
    ConstantCoefficient coeff_temperature_in(InflowTemperature);
    b0->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(coeff_temperature_in, coeff_rA0, config.sigma, config.kappa), ess_bdr_0);
    b0->AddBdrFaceIntegrator(new BoundaryFlowIntegrator(coeff_temperature_in, coeff_rV, -1.), ess_bdr_0);
    b0->Assemble();
    B0 = b0->ParallelAssemble();

    if (b1) delete b1;
    if (B1) delete B1;
    b1 = new ParLinearForm(&fespace_L2);
    ConstantCoefficient coeff_salinity_in(InflowSalinity);
    b1->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(coeff_salinity_in, coeff_rA1, config.sigma, config.kappa), ess_bdr_1);
    b1->AddBdrFaceIntegrator(new BoundaryFlowIntegrator(coeff_salinity_in, coeff_rV, -1.), ess_bdr_1);
    b1->Assemble();
    B1 = b1->ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Recover actual information
    ParGridFunction temperature(&fespace_L2);
    ParGridFunction salinity(&fespace_L2);

    temperature.Distribute(X.GetBlock(0));
    salinity.Distribute(X.GetBlock(1));

    //Create the auxiliar grid functions
    ParGridFunction inverse_permeability(&fespace_L2);
    ParGridFunction expansivity_temperature(&fespace_L2);
    ParGridFunction expansivity_salinity(&fespace_L2);
    ParGridFunction temperature_dr(&fespace_L2);
    ParGridFunction salinity_dr(&fespace_L2);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < temperature.Size(); ii++){
        inverse_permeability(ii) = InversePermeability(temperature(ii), salinity(ii));
        expansivity_temperature(ii) = ExpansivityTemperature(temperature(ii), salinity(ii));
        expansivity_salinity(ii) = ExpansivitySalinity(temperature(ii), salinity(ii));
    }

    temperature.GetDerivative(1, 0, temperature_dr);
    salinity.GetDerivative(1, 0, salinity_dr);

    //Set the grid coefficients
    GridFunctionCoefficient coeff_gamma(&inverse_permeability);

    GridFunctionCoefficient coeff_expansivity_temperature(&expansivity_temperature);
    GridFunctionCoefficient coeff_expansivity_salinity(&expansivity_salinity);
    GridFunctionCoefficient coeff_temperature_dr(&temperature_dr);
    GridFunctionCoefficient coeff_salinity_dr(&salinity_dr);

    ProductCoefficient coeff_neg_gamma(-1., coeff_gamma);
    ScalarVectorProductCoefficient coeff_neg_r_inv_gamma_hat(coeff_neg_gamma, coeff_r_inv_hat);

    ProductCoefficient coeff_buoyancy_temperature(coeff_expansivity_temperature, coeff_temperature_dr);
    ProductCoefficient coeff_buoyancy_salinity(coeff_expansivity_salinity, coeff_salinity_dr);
    ProductCoefficient coeff_r_buoyancy_temperature(coeff_r, coeff_buoyancy_temperature);
    ProductCoefficient coeff_r_buoyancy_salinity(coeff_r, coeff_buoyancy_salinity);
    
    //Define non-constant bilinear forms of the system
    if(A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_gamma));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_r_inv_gamma_hat));
    a11.Assemble();
    a11.Finalize();
    A11 = a11.ParallelAssemble();
    A11_e = A11->EliminateRowsCols(ess_tdof_1);

    //Define the non-constant RHS
    ParLinearForm b1(&fespace_H1);
    //b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_temperature));
    //b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_buoyancy_salinity));
    b1.Assemble();
    b1.ParallelAssemble(B1);
    A10_e->Mult(Vorticity, B1, -1., 1.);
    EliminateBC(*A11, *A11_e, ess_tdof_1, Stream, B1);
}
