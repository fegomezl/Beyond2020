#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    transport_oper->SetParameters(X, *rV);
    ode_solver->Step(X, t, dt);

    flow_oper->SetParameters(X);
    flow_oper->Solve(Z, *V, *rV);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        theta->Distribute(&(X.GetBlock(0)));
        phi->Distribute(&(X.GetBlock(1)));
        w->Distribute(&(Z.GetBlock(0)));
        psi->Distribute(&(Z.GetBlock(1)));
        v->Distribute(V);
        rv->Distribute(rV);

        //Calculate phases
        for (int ii = 0; ii < phase->Size(); ii++){
            double T_f = config.T_f + T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*EpsilonInv*((*theta)(ii) - T_f)));
        }

        //Normalize stream
        if (config.rescale){
            double psi_local_max = psi->Max(), psi_max;
            double psi_local_min = psi->Min(), psi_min;
            MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            for (int ii = 0; ii < psi->Size(); ii++)
                    (*psi)(ii) = ((*psi)(ii)-psi_min)/(psi_max-psi_min);
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
    theta.SetFromTrueDofs(X.GetBlock(0));
    phi.SetFromTrueDofs(X.GetBlock(1));

    theta.GetDerivative(1, 0, theta_dr);
    phi.GetDerivative(1, 0, phi_dr);

    //Calculate eta and buoyancy coefficients
    for (int ii = 0; ii < phase.Size(); ii++){
        double T = theta(ii);
        double S = phi(ii);
        double P = 0.5*(1 + tanh(5*EpsilonInv*(theta(ii) - config.T_f - T_fun(phi(ii))))); 

        eta(ii) = Epsilon + pow(1-P, 2)/(pow(P, 3) + Epsilon);

        theta(ii) = delta_rho_t_fun(T, S);
        phi(ii) = delta_rho_p_fun(T, S);
    }

    //Properties coefficients
    GridFunctionCoefficient Eta(&eta);

    GridFunctionCoefficient Theta_dr(&theta_dr);
    GridFunctionCoefficient k_t(&theta);
    ProductCoefficient k_Theta_dr(k_t, Theta_dr);

    GridFunctionCoefficient Phi_dr(&phi_dr);
    GridFunctionCoefficient k_p(&phi);
    ProductCoefficient k_Phi_dr(k_p, Phi_dr);

    //Rotational coupled coefficients
    ScalarVectorProductCoefficient Eta_r_inv_hat(Eta, r_inv_hat);
    ProductCoefficient k_r_Theta_dr(coeff_r, k_Theta_dr);
    ProductCoefficient k_r_Phi_dr(coeff_r, k_Phi_dr);

    //Define non-constant bilinear forms of the system
    if (D) delete D;
    if (D_e) delete D_e;
    ParBilinearForm d(&fespace);
    d.AddDomainIntegrator(new DiffusionIntegrator(Eta));
    d.AddDomainIntegrator(new ConvectionIntegrator(Eta_r_inv_hat));
    d.Assemble();
    d.Finalize();
    D = d.ParallelAssemble();
    *D *= -1;
    D_e = D->EliminateRowsCols(ess_tdof_psi);

    //Define the non-constant RHS
    ParLinearForm f(&fespace);
    f.AddDomainIntegrator(new DomainLFIntegrator(k_r_Theta_dr));
    f.AddDomainIntegrator(new DomainLFIntegrator(k_r_Phi_dr));
    f.Assemble();
    f.ParallelAssemble(F);
    Ct_e->Mult(W, F, -1., 1.);
    EliminateBC(*D, *D_e, ess_tdof_psi, Psi, F);
}
