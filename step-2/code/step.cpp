#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t >= (config.t_final + config.t_init) - 1e-8*config.dt_init);
    dt = min(dt, (config.t_final + config.t_init) - t);

    //Update boundary conditions
    boundary.SetTime(t);
    x->ProjectBdrCoefficient(boundary, ess_bdr);
    x->GetTrueDofs(X);

    //Perform the time_step
    oper->SetParameters(X, ess_bdr, *x);
    ode_solver->Step(X, t, dt);

    //Output from the solution
    if (last || (iteration % config.vis_steps) == 0){
        //Calculate convergence
        x->SetFromTrueDofs(X);
        actual_error = x->ComputeL2Error(boundary);
        total_error += actual_error;
        iterations_error += 1;

        //Graph
        paraview_out->SetCycle(iteration);
        paraview_out->SetTime(t - config.t_init);
        paraview_out->Save();
    }
    
    //Print the system state 
    double percentage = 100*(t - config.t_init)/config.t_final;
    string progress = to_string((int)percentage)+"%";
    if (config.master){    
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << t - config.t_init << setw(12)
             << progress << setw(12)
             << actual_error << "\r";
        cout.flush();
    }
}

void Conduction_Operator::SetParameters(const Vector &X, Array<int> ess_bdr, ParGridFunction &x){
    //Read the solution X
    ParGridFunction aux(&fespace);
    aux.SetFromTrueDofs(X);

    //Create the K coefficient
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > T_f)
            aux(ii) *= alpha_l;
        else  
            aux(ii) *= alpha_s;
    }
    GridFunctionCoefficient coeff(&aux);

    //Construct M
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble(0);
    m->FormLinearSystem(ess_tdof_list, x, *f, M, Z, F, 1);
    M_solver.SetOperator(M);

    //Create the new K
    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff));
    k->Assemble(0);
    m->FormLinearSystem(ess_tdof_list, x, *f, K, Z, F, 1);
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    K.Mult(X,z);
    z.Neg();
    M_solver.Mult(z, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //Solve M(dX_dt) = -K(X + dt*dX_dt)] for dX_dt
    if (T) delete T;
    T = Add(1., M, dt, K);
    T_solver.SetOperator(*T);

    K.Mult(X, z);
    z.Neg();
    T_solver.Mult(z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &b, 
                             int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    if (T) delete T;
    T = Add(1., M, scaled_dt, K);
    T_solver.SetOperator(*T);
    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &b, Vector &X, 
                             double tol){
    //Solve the system Ax = z -> (M - gamma*K)x = Mb
    M.Mult(b,z);
    T_solver.Mult(z,X);
    return 0;
}
