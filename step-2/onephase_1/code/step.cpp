#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    boundary.SetTime(t);
    x->ProjectBdrCoefficient(boundary, ess_bdr);
    x->GetTrueDofs(*X);

    oper->SetParameters(*X);
    ode_solver->Step(*X, t, dt);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Graph
        x->SetFromTrueDofs(*X);
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
    aux.SetFromTrueDofs(X);

    //Create the K coefficient
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > T_f)
            aux(ii) = alpha_l;
        else
            aux(ii) = alpha_s;
    }
    alpha.SetGridFunction(&aux);
    r_alpha.SetACoef(r); r_alpha.SetBCoef(alpha);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(r_alpha));
    k->Assemble();
    k->Finalize();
}

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    ParGridFunction tmp_X(&fespace);
    ParLinearForm tmp_z(&fespace);
    tmp_X.SetFromTrueDofs(X);
    k->Mult(tmp_X,tmp_z);
    tmp_z.Neg();

    OperatorHandle A;
    HypreParVector B(&fespace);
    ParGridFunction tmp_dX_dt(&fespace);
    tmp_dX_dt = 0.0;

    m->FormLinearSystem(ess_tdof_list, tmp_dX_dt, tmp_z, A, dX_dt, B);
    M_solver.Mult(B, dX_dt);
}

void Conduction_Operator::ImplicitSolve(const double dt, const Vector &X, Vector &dX_dt){
    //Solve M(dX_dt) = -K(X + dt*dX_dt)] for dX_dt
    if(t) delete t;
    t = new ParBilinearForm(&fespace);
    ParGridFunction tmp_X(&fespace);
    ParLinearForm tmp_z(&fespace);
    tmp_X.SetFromTrueDofs(X);
    t->AddDomainIntegrator(new MassIntegrator(r));
    dt_r_alpha.SetAConst(dt); dt_r_alpha.SetBCoef(r_alpha);
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_r_alpha));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T_op);
    T_solver.SetOperator(T_op);

    k->Mult(tmp_X, tmp_z);
    tmp_z.Neg();

    OperatorHandle A;
    Vector B;
    ParGridFunction tmp_dX_dt(&fespace);
    tmp_dX_dt = 0.0;

    t->FormLinearSystem(ess_tdof_list, tmp_dX_dt, tmp_z, A, dX_dt, B);

    T_solver.Mult(B, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &b,
                             int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    if(t) delete t;
    t = new ParBilinearForm(&fespace);
    SUN_tmp_X = ParGridFunction(&fespace);
    SUN_tmp_z = ParLinearForm(&fespace);
    SUN_tmp_X.SetFromTrueDofs(X);
    t->AddDomainIntegrator(new MassIntegrator(r));
    dt_r_alpha.SetAConst(scaled_dt); dt_r_alpha.SetBCoef(r_alpha);
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_r_alpha));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T_op);
    T_solver.SetOperator(T_op);

    OperatorHandle A;
    Vector SUN_X;

    t->FormLinearSystem(ess_tdof_list, SUN_tmp_X, SUN_tmp_z, A, SUN_X, SUN_B);

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
