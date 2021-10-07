#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    ode_solver->Step(*X, t, dt);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Update information
        x->SetFromTrueDofs(*X);

        //Convergence analysis
        exact.SetTime(t);
        actual_error = x->ComputeL2Error(exact, irs)/ComputeGlobalLpNorm(2, exact, *pmesh, irs); 
        total_error += actual_error;

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
             << progress << setw(9)
             << actual_error << "\r";
        cout.flush();
    }
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X) = F for dX_dt
    
    Z = 0.;
    dX_dt = 0.;
    
    K_0->Mult(-1., X, 1., Z);
    EliminateBC(*M, *M_e, ess_tdof_list, dX_dt, Z);

    M_solver.Mult(Z, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //From  M(dX_dt) + K(X) = F
    //Solve M(dX_dt) + K(X + dt*dX_dt) = F for dX_dt 
    
    Z = 0.;
    dX_dt = 0.;

    if (T) delete T;
    if (T_e) delete T_e;
    T = Add(1., *M_0, dt, *K_0);
    T_e = T->EliminateRowsCols(ess_tdof_list); 
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

    K_0->Mult(-1., X, 1., Z);
    EliminateBC(*T, *T_e, ess_tdof_list, dX_dt, Z);

    T_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
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

int Conduction_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new

    Z = 0.;
    X_new = X;

    M_0->Mult(X, Z);
    EliminateBC(*T, *T_e, ess_tdof_list, X_new, Z);

    T_solver.Mult(Z, X_new);
    return 0;
}
