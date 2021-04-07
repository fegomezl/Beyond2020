#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t + dt >= config.t_final - dt/2.);

    //Perform the time_step
    oper->SetParameters(X);
    ode_solver->Step(X, t, dt);

    //Print the system state, if neccesary
    if (last || (iteration % config.vis_steps) == 0){
        double percentage = 100*(t - config.t_init)/(config.t_final - config.t_init);
        string progress = to_string((int)percentage)+"%";
        if (config.master) 
            cout << left << setw(8)
                 << iteration << setw(8)
                 << t << setw(8)
                 << progress << "\n";
        x->SetFromTrueDofs(X);
        paraview_out->SetCycle(iteration);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }
}

void Conduction_Operator::SetParameters(const Vector &X){
    //Initialize the bilinear forms
    delete k;

    //Define K coefficient
    ParGridFunction x_k(fespace);
    x_k.SetFromTrueDofs(X);
    for (int ii = 0; ii < x_k.Size(); ii++){
        if (x_k(ii) < T_f) x_k(ii) *= alpha_s;
        else x_k(ii) *= alpha_l;
    }
    GridFunctionCoefficient coeff_k(&x_k);

    //Create the new K
    k = new ParBilinearForm(fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_k));
    k->Assemble(0);
    k->FormSystemMatrix(ess_tdof_list, K);

    //Renovate T
    delete T;
    T = NULL;

    //Update the solvers
    M_solver.SetOperator(M);
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
