#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(*X, *r_Velocity);
    ode_solver->Step(*X, t, dt);

    flow_oper->SetParameters(*X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        enthalphy->Distribute(X);
        vorticity->Distribute(Y.GetBlock(0));
        stream->Distribute(Y.GetBlock(1));
        velocity->Distribute(Velocity);
        r_velocity->Distribute(r_Velocity);

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

void Transport_Operator::SetParameters(const Vector &X, const Vector &r_Velocity){
    //Create the auxiliar grid functions
    ParGridFunction enthalphy(&fespace_L2);
    ParGridFunction r_velocity(&fespace_ND);

    enthalphy.Distribute(X);
    r_velocity.Distribute(r_Velocity);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < enthalphy.Size(); ii++){
        //Map option
        double H = enthalphy(ii);
        enthalphy(ii) = config.a_s*(0.5*(1 - tanh(5*JumpScale*(H))))
                + config.a_l*(0.5*(1 + tanh(5*JumpScale*(H))));
        //aux(ii) = config.a_s*(0.5*(1 - tanh(5*JumpScale*(H+1))))
        //        + config.a_l*(0.5*(1 + tanh(5*JumpScale*(H-1))));
    }

    //Set the grid coefficients
    enthalphy.ExchangeFaceNbrData();
    GridFunctionCoefficient coeff_A(&enthalphy);
    ProductCoefficient coeff_rA(coeff_r, coeff_A);
    VectorGridFunctionCoefficient coeff_rV(&r_velocity);
    
    config.sigma = 1.;
    config.kappa = 0.;
    config.eta   = 0.;

    //Create RHS
    ConstantCoefficient zero(0.);
    if (b) delete b;
    if (B) delete B;
    b = new ParLinearForm(&fespace_L2);
    b->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(zero, coeff_rA, config.sigma, config.kappa), ess_bdr);
    b->Assemble();
    B = b->ParallelAssemble();

    //Create transport matrix
    if (k) delete k;
    if (K) delete K;
    k = new ParBilinearForm(&fespace_L2);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rA, config.sigma, config.kappa));
    k->AddInteriorFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.));
    k->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA, config.sigma, config.kappa), ess_bdr);
    k->AddBdrFaceIntegrator(new NonconservativeDGTraceIntegrator(coeff_rV, 1.), ess_bdr);

    //Dont know how it works
    //k->AddInteriorFaceIntegrator(new DGDiffusionBR2Integrator(&fespace, config.eta));
    //k->AddBdrFaceIntegrator(new DGDiffusionBR2Integrator(&fespace, config.eta), ess_bdr);

    k->Assemble();
    k->Finalize();
    K = k->ParallelAssemble();    
}

void Flow_Operator::SetParameters(const Vector &X){

    //Recover actual information
    ParGridFunction enthalphy(&fespace_L2), density(&fespace_L2);
    ParGridFunction density_dr(&fespace_L2);
    ParGridFunction permeability(&fespace_L2);

    enthalphy.Distribute(X);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < permeability.Size(); ii++){
        double P = Phase(enthalphy(ii));

        permeability(ii) = Epsilon + pow(1-P, 2)/(pow(P, 3) + Epsilon);
        density(ii) = Density(enthalphy(ii));
    }

    density.GetDerivative(1, 0, density_dr);

    //Properties coefficients
    GridFunctionCoefficient coeff_permeability(&permeability);
    ProductCoefficient coeff_neg_permeability(-1., coeff_permeability);

    GridFunctionCoefficient coeff_density_dr(&density_dr);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient coeff_neg_r_inv_permeability_hat(coeff_neg_permeability, coeff_r_inv_hat);
    ProductCoefficient coeff_r_density_dr(coeff_r, coeff_density_dr);

    //Define non-constant bilinear forms of the system
    if(A11) delete A11;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_neg_permeability));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_neg_r_inv_permeability_hat));
    a11.Assemble();
    a11.Finalize();
    A11 = a11.ParallelAssemble();
    A11_e = A11->EliminateRowsCols(ess_tdof_1);

    //Define the non-constant RHS
    ParLinearForm b1(&fespace_H1);
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_r_density_dr));
    b1.Assemble();
    b1.ParallelAssemble(B1);
    A10_e->Mult(Vorticity, B1, -1., 1.);
    EliminateBC(*A11, *A11_e, ess_tdof_1, Stream, B1);
}
