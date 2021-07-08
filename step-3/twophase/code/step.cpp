#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    oper_T->SetParameters(*X_T);
    ode_solver->Step(*X_T, t, dt);

    //Solve the flow problem
    flow_oper->Solve(config, X_Psi, x_psi, X_T);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Graph
        x_T->SetFromTrueDofs(*X_T);
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

void Conduction_Operator::SetParameters(const Vector &X){
    //Create the auxiliar grid functions

    Vector Aux(X);
    if (config.T_f != 0)
        Aux -= config.T_f;
    aux.SetFromTrueDofs(Aux);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
        } else {
            aux_C(ii) = config.c_s;
            aux_K(ii) = config.k_s;
        }

        aux(ii) = config.L*config.invDeltaT*exp(-M_PI*pow(config.invDeltaT*aux(ii), 2));
    }

    //Set the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_L(&aux);

    SumCoefficient coeff_CL(coeff_C, coeff_L);

    coeff_rK.SetBCoef(coeff_K); dt_coeff_rK.SetBCoef(coeff_rK); 
    coeff_rCL.SetBCoef(coeff_CL);
    coeff_rCLV.SetACoef(coeff_CL);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

    //Create corresponding bilinear forms
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    m->Assemble();
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k->Assemble();
    k->Finalize();
}

void Conduction_Operator::UpdateVelocity(const HypreParVector &Psi){
    psi.SetFromTrueDofs(Psi);
    gradpsi.SetGridFunction(&psi);
    rV.SetBCoef(gradpsi);
    coeff_rCLV.SetBCoef(rV);
}
