#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t + dt >= config.t_final - dt/2.);

    //Perform the time_step
    oper->SetParameters(X);
    ode_solver->Step(X, t, dt);

    //Print the system state 
    double percentage = 100*t/config.t_final;
    string progress = to_string((int)percentage)+"%";
    if (config.master){
        cout << left << setw(8)
             << iteration << setw(8)
             << t << setw(8)
             << progress << "\r";
        cout.flush();
    }
    if (last || (iteration % config.vis_steps) == 0){
        x->SetFromTrueDofs(X);
        paraview_out->SetCycle(iteration);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }
}

void Conduction_Operator::SetParameters(const Vector &X){
    //Read the solution x
    ParGridFunction x(fespace);
    x.SetFromTrueDofs(X);

    //Create the K coefficient
    for (int ii = 0; ii < x.Size(); ii++){
        if (x(ii) > T_f) 
            x(ii) *= alpha_l;
        else  
            x(ii) *= alpha_s;
    }
    GridFunctionCoefficient coeff(&x);

    //Create the new K
    delete k;
    k = new ParBilinearForm(fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff));
    k->Assemble(0);
    k->FormSystemMatrix(ess_tdof_list, K);

    //Renovate T
    delete T;
    T = NULL;
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    K.Mult(X,z);
    z.Neg();
    M_solver.Mult(z, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //Solve M(dX_dt) = -K(X + dt*dX_dt)] for dX_dt
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
