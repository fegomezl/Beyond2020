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
        vis_impressions++;

        //Update information
        temperature->Distribute(X.GetBlock(0));
        salinity->Distribute(X.GetBlock(1));
        vorticity->Distribute(Y.GetBlock(0));
        stream->Distribute(Y.GetBlock(1));
        velocity->Distribute(Velocity);
        rvelocity->Distribute(rVelocity);

        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = config.T_f + T_fun((*salinity)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*EpsilonInv*((*temperature)(ii) - T_f)));
        }

        //Normalize stream
        if (config.rescale){
            double psi_local_max = stream->Max(), psi_max;
            double psi_local_min = stream->Min(), psi_min;
            MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            for (int ii = 0; ii < stream->Size(); ii++)
                    (*stream)(ii) = ((*stream)(ii)-psi_min)/(psi_max-psi_min);
        }

        //Graph
        paraview_out->SetCycle(vis_impressions);
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
        double DT = temperature(ii) - config.T_f - T_fun(salinity(ii));
        double P = 0.5*(1 + tanh(5*EpsilonInv*DT));

        aux_C(ii) = config.c_s + (config.c_l-config.c_s)*P;
        aux_K(ii) = config.k_s + (config.k_l-config.k_s)*P;
        aux_D(ii) = config.d_s + (config.d_l-config.d_s)*P;

        temperature(ii) = DT;
        phase(ii) = P;
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);

    //Construct latent heat term
    GradientGridFunctionCoefficient coeff_dT(&temperature);
    GradientGridFunctionCoefficient coeff_dH(&phase);
    
    InnerProductCoefficient coeff_dHdT(coeff_dH, coeff_dT);
    InnerProductCoefficient coeff_dT_2(coeff_dT, coeff_dT);

    SumCoefficient coeff_dT_2e(Epsilon, coeff_dT_2);
    
    PowerCoefficient coeff_inv_dT_2(coeff_dT_2e, -1.);
    ProductCoefficient coeff_DeltaT(coeff_dHdT, coeff_inv_dT_2);
    ProductCoefficient coeff_LDeltaT(0.5*(config.L_s+config.L_l), coeff_DeltaT);

    SumCoefficient coeff_CL(coeff_C, coeff_LDeltaT);

    //Construct final coefficients
    ProductCoefficient coeff_rCL(coeff_r, coeff_CL);
    ProductCoefficient coeff_rK(coeff_r, coeff_K);
    ProductCoefficient coeff_rD(coeff_r, coeff_D);

    VectorGridFunctionCoefficient coeff_rV(&rvelocity);
    ScalarVectorProductCoefficient coeff_rCLV(coeff_CL, coeff_rV);

    //Create corresponding bilinear forms
    if (M0) delete M0;
    if (M0_e) delete M0_e;
    if (M0_o) delete M0_o;
    ParBilinearForm m0(&fespace_H1);
    m0.AddDomainIntegrator(new MassIntegrator(coeff_rCL));
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
    k0.AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k0.AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k0.Assemble();
    k0.Finalize();
    K0 = k0.ParallelAssemble();    

    if (K1) delete K1;
    ParBilinearForm k1(&fespace_H1);
    k1.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k1.AddDomainIntegrator(new ConvectionIntegrator(coeff_rV));
    k1.Assemble();
    k1.Finalize();
    K1 = k1.ParallelAssemble();
}

void Flow_Operator::SetParameters(const BlockVector &X){
    //Update information
    temperature.SetFromTrueDofs(X.GetBlock(0));
    salinity.SetFromTrueDofs(X.GetBlock(1));

    temperature.GetDerivative(1, 0, temperature_dr);
    salinity.GetDerivative(1, 0, salinity_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < phase.Size(); ii++){
        double T = temperature(ii);
        double S = salinity(ii);
        double P = 0.5*(1 + tanh(5*EpsilonInv*(T - config.T_f - T_fun(S)))); 

        impermeability(ii) = Epsilon + pow(1-P, 2)/(pow(P, 3) + Epsilon);

        temperature(ii) = delta_rho_t_fun(T, S);
        salinity(ii) = delta_rho_p_fun(T, S);
    }

    //Properties coefficients
    GridFunctionCoefficient coeff_impermeability(&impermeability);

    GridFunctionCoefficient coeff_temperature_dr(&temperature_dr);
    GridFunctionCoefficient coeff_k_t(&temperature);
    ProductCoefficient coeff_k_temperature_dr(coeff_k_t, coeff_temperature_dr);

    GridFunctionCoefficient coeff_salinity_dr(&salinity_dr);
    GridFunctionCoefficient coeff_k_p(&salinity);
    ProductCoefficient coeff_k_salinity_dr(coeff_k_p, coeff_salinity_dr);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient coeff_impermeability_r_inv_hat(coeff_impermeability, coeff_r_inv_hat);
    ProductCoefficient coeff_k_r_temperature_dr(coeff_r, coeff_k_temperature_dr);
    ProductCoefficient coeff_k_r_salinity_dr(coeff_r, coeff_k_salinity_dr);

    //Define non-constant bilinear forms of the system
    if (A11) delete A11;
    if (A11_e) delete A11_e;
    ParBilinearForm a11(&fespace_H1);
    a11.AddDomainIntegrator(new DiffusionIntegrator(coeff_impermeability));
    a11.AddDomainIntegrator(new ConvectionIntegrator(coeff_impermeability_r_inv_hat));
    a11.Assemble();
    a11.Finalize();
    A11 = a11.ParallelAssemble();
    *A11 *= -1;
    A11_e = A11->EliminateRowsCols(ess_tdof_1);

    //Define the non-constant RHS
    if (B1) delete B1;
    ParLinearForm b1(&fespace_H1);
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_k_r_temperature_dr));
    b1.AddDomainIntegrator(new DomainLFIntegrator(coeff_k_r_salinity_dr));
    b1.Assemble();
    B1 = b1.ParallelAssemble();
    A10_e->Mult(Vorticity, *B1, -1., 1.);
    EliminateBC(*A11, *A11_e, ess_tdof_1, Stream, *B1);
}
