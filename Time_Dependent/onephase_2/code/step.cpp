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
    dX_dt = 0.;
    
    HypreParVector Z(&fespace);
    K0.Mult(-1., X, 1., Z);

    EliminateBC(M, *Me, ess_tdof_list, dX_dt, Z);

    M_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K
    if (Te) delete Te;
    T = *Add(1., M0, scaled_dt, K0);
    Te = m->ParallelEliminateTDofs(ess_tdof_list, T);
    
    T_prec.SetOperator(T);
    T_solver.SetOperator(T);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = F
    //Solve M(X_new - X) + dt*K(X_new) = dt*F for X_new
    X_new = X;
    HypreParVector Z(&fespace);
    M0.Mult(X, Z);

    EliminateBC(T, *Te, ess_tdof_list, X_new, Z);

    T_solver.Mult(Z, X_new);
    return 0;
}
