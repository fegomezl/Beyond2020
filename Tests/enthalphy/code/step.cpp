#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    oper->SetParameters(*X);
    ode_solver->Step(*X, t, dt);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        x->SetFromTrueDofs(*X);

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

void Transport_Operator::SetParameters(const Vector &X){
    //Create the auxiliar grid functions
    aux.SetFromTrueDofs(X);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        //Map option
        double H = aux(ii);
        aux(ii) = config.a_s*(0.5*(1 - tanh(5*JumpScale*(H+1))))
                + config.a_l*(0.5*(1 + tanh(5*JumpScale*(H-1))));

        //Conditional option
        /*
        if (aux(ii) >= 1)
            aux(ii) = config.a_l;
        else if (aux(ii) <= -1)
            aux(ii) = config.a_s;
        else
            aux(ii) = 0.;
        */
    }

    //Set the grid coefficients
    GridFunctionCoefficient coeff_A(&aux);
    coeff_rA.SetBCoef(coeff_A);

    //Create transport matrix
    if (k) delete k;
    if (K_0) delete K_0;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA));
    k->Assemble();
    k->Finalize();
    K_0 = k->ParallelAssemble();    
}

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt
    
    Z = 0.;
    dX_dt = 0.;
    
    K_0->Mult(-1., X, 1., Z);
    EliminateBC(*M, *M_e, ess_tdof_list, dX_dt, Z);

    M_solver.Mult(Z, dX_dt);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K

    if (T) delete T;
    if (T_e) delete T_e;
    T = Add(1., *M_0, scaled_dt, *K_0);
    T_e = T->EliminateRowsCols(ess_tdof_list);
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    Z = 0.;
    X_new = X;

    M_0->Mult(X, Z);
    EliminateBC(*T, *T_e, ess_tdof_list, X_new, Z);

    T_solver.Mult(Z, X_new);
    return 0;
}
