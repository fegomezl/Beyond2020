#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t + dt >= config.t_final - dt/2.);

    //Perform the time_step
    ode_solver->Step(X, t, dt);

    //Print the system state, if neccesary
    if (last || (iteration % config.vis_steps) == 0){
        if (config.master) 
            cout << "step " << iteration 
                 << ", t = " << t << ".\n";
        x->SetFromTrueDofs(X);
        paraview_out->SetCycle(iteration);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }

    //Update the operator
    oper->SetParameters(X);
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Compute dX_dt = M^-1*-K(X) for dX_dt
    K.Mult(X,z);
    z.Neg();
    M_solver.Mult(z, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //Solve the equation dX_dt = M^-1*[-K(X + dt*dX_dt)] for dX_dt
    if (!T){
        T = Add(1.0, M, dt, K);
        current_dt = dt;
        T_solver.SetOperator(*T);
    }
    MFEM_VERIFY(dt == current_dt, ""); //?
    K.Mult(X, z);
    z.Neg();
    T_solver.Mult(z, dX_dt);

}
