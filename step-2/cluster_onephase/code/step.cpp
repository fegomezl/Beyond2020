#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
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

        //Calculate convergence
        initial_f.SetTime(t);
        x->SetFromTrueDofs(*X);
        actual_error = x->ComputeL2Error(initial_f);
        actual_error /= Zmax*Rmax;
        total_error += actual_error;

        if (config.master && config.print){
            ofstream output;
            string fname = "results/local_" + to_string(config.ode_solver_type) + "_" + to_string(config.order) + "_" + to_string(config.refinements) + ".txt";
            output.open(fname, ios::app);
            output.precision(4);
            output << left << setw(12)
                   << t << setw(12)
                   << actual_error << "\n";
            output.close();
        }
    }
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    ParGridFunction x(&fespace);
    ParLinearForm z(&fespace);
    x.SetFromTrueDofs(X);
    k->Mult(x, z);
    z.Neg();

    OperatorHandle A;
    HypreParVector Z(&fespace);
    ParGridFunction dx_dt(&fespace);
    dx_dt = 0.;

    m->FormLinearSystem(ess_tdof_list, dx_dt, z, A, dX_dt, Z);
    M_solver.Mult(Z, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //Solve M(dX_dt) = -K(X + dt*dX_dt)] for dX_dt
    if (t) delete t;
    t = new ParBilinearForm(&fespace);
    t->AddDomainIntegrator(new MassIntegrator(r));
    dt_r_alpha.SetAConst(dt);
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_r_alpha));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T);
    T_solver.SetOperator(T);

    ParGridFunction x(&fespace);
    ParLinearForm z(&fespace);
    x.SetFromTrueDofs(X);
    k->Mult(x, z);
    z.Neg();

    OperatorHandle A;
    HypreParVector Z(&fespace);
    ParGridFunction dx_dt(&fespace);
    dx_dt = 0.;

    t->FormLinearSystem(ess_tdof_list, dx_dt, z, A, dX_dt, Z);
    T_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    if (t) delete t;
    t = new ParBilinearForm(&fespace);
    t->AddDomainIntegrator(new MassIntegrator(r));
    dt_r_alpha.SetAConst(scaled_dt);
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_r_alpha));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T);
    T_solver.SetOperator(T);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &B, Vector &X, double tol){
    //Solve the system Ax = z -> (M - gamma*K)x = Mb
    HypreParVector Z(&fespace);
    M.Mult(B,Z);
    T_solver.Mult(Z,X);
    return 0;
}
