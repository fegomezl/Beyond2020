#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    cond_oper->SetParameters(X);
    ode_solver->Step(X, t, dt);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        theta->Distribute(&(X.GetBlock(0)));
        phi->Distribute(&(X.GetBlock(1)));

        //Calculate phases
        double T_f;
        for (int ii = 0; ii < phase->Size(); ii++){
            T_f = T_fun((*phi)(ii));
            (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
        }

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

void Conduction_Operator::SetParameters(const BlockVector &X){
    //Recover actual information
    theta.Distribute(&(X.GetBlock(0)));
    phi.Distribute(&(X.GetBlock(1)));

    //Associate the values of each auxiliar function
    double DT = 0.;
    for (int ii = 0; ii < phase.Size(); ii++){
        DT = theta(ii) - T_fun(phi(ii)); 
        theta(ii) = DT;
        phase(ii) = 0.5*(1 + tanh(5*config.invDeltaT*DT));
    }

    aux_C.Set(config.c_l-config.c_s, phase); aux_C += config.c_s;
    aux_K.Set(config.k_l-config.k_s, phase); aux_K += config.k_s;
    aux_D.Set(config.d_l-config.d_s, phase); aux_D += config.d_s;
    aux_L.Set(config.L_l-config.L_s, phase); aux_L += config.L_s;

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_D(&aux_D);
    GridFunctionCoefficient coeff_L(&aux_L);

    //Construct latent heat term
    GradientGridFunctionCoefficient dT(&theta);
    GradientGridFunctionCoefficient dH(&phase);
    
    dHdT.SetACoef(dH);  dT_2.SetACoef(dT);
    dHdT.SetBCoef(dT);  dT_2.SetBCoef(dT);

    SumCoefficient dT_2e(config.EpsilonT, dT_2);
    
    PowerCoefficient inv_dT_2(dT_2e, -1.);
    ProductCoefficient DeltaT(dHdT, inv_dT_2);
    ProductCoefficient LDeltaT(coeff_L, DeltaT);

    SumCoefficient coeff_CL(coeff_C, LDeltaT);

    //Construct final coefficients
    coeff_rC.SetBCoef(coeff_CL);
    coeff_rK.SetBCoef(coeff_K); 
    coeff_rD.SetBCoef(coeff_D); 

    //Create corresponding bilinear forms
    ParBilinearForm m_theta(&fespace);
    HypreParMatrix *M_theta;
    m_theta.AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m_theta.Assemble();
    m_theta.Finalize();
    M_theta = m_theta.ParallelAssemble();

    ParBilinearForm m_phi(&fespace);
    HypreParMatrix *M_phi;
    m_phi.AddDomainIntegrator(new MassIntegrator(coeff_r));
    m_phi.Assemble();
    m_phi.Finalize();
    M_phi = m_phi.ParallelAssemble();

    ParBilinearForm k_theta(&fespace);
    HypreParMatrix *K_theta;
    k_theta.AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k_theta.Assemble();
    k_theta.Finalize();
    K_theta = k_theta.ParallelAssemble();

    ParBilinearForm k_phi(&fespace);
    HypreParMatrix *K_phi;
    k_phi.AddDomainIntegrator(new DiffusionIntegrator(coeff_rD));
    k_phi.Assemble();
    k_phi.Finalize();
    K_phi = k_phi.ParallelAssemble();

    //Set solver objects 
    delete M_0;
    Array2D<HypreParMatrix*> MBlocks(2,2);
    MBlocks = NULL;
    MBlocks(0, 0) = M_theta;
    MBlocks(1, 1) = M_phi;
    M_0 = HypreParMatrixFromBlocks(MBlocks);

    delete K_0;
    Array2D<HypreParMatrix*> KBlocks(2,2);
    KBlocks = NULL;
    KBlocks(0, 0) = K_theta;
    KBlocks(1, 1) = K_phi;
    K_0 = HypreParMatrixFromBlocks(KBlocks);

    delete M;
    HypreParMatrix M_aux(*M_0);
    M_aux.EliminateRowsCols(ess_tdof_list);
    M = new SuperLURowLocMatrix(M_aux);

    //if (M_solver)
    //    M_solver->DismantleGrid();
    delete M_solver;    
    M_solver = new SuperLUSolver(MPI_COMM_WORLD);
    M_solver->SetPrintStatistics(false);
    M_solver->SetSymmetricPattern(true);
    M_solver->SetColumnPermutation(superlu::PARMETIS);
    M_solver->SetIterativeRefine(superlu::SLU_DOUBLE);
    M_solver->SetOperator(*M);

    //Delete used memory
    delete M_theta;
    delete M_phi;
    delete K_theta;
    delete K_phi;
}
