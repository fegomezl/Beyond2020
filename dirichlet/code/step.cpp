#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    boundary.SetTime(t);
    x->ProjectBdrCoefficient(boundary, ess_bdr);

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
    ParGridFunction aux(&fespace);
    aux.SetFromTrueDofs(X);

    //Create the K coefficient
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > T_f)
            aux(ii) = alpha_l;
        else
            aux(ii) = alpha_s;
    }
    GridFunctionCoefficient alpha(&aux);
    ProductCoefficient r_alpha(alpha, r);
    if(r_alpha_dt) delete r_alpha_dt;
    r_alpha_dt = new ProductCoefficient(1.,r_alpha);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(r_alpha));
    k->Assemble();
    k->Finalize();
    //k->FormSystemMatrix(ess_tdof_list, K);
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
    //if (T) delete T;
    ParBilinearForm t(&fespace);
    ParGridFunction tmp_X(&fespace);
    ParLinearForm tmp_z(&fespace);
    HypreParMatrix T_op;
    tmp_X.SetFromTrueDofs(X);
    t.AddDomainIntegrator(new MassIntegrator(r));
    r_alpha_dt->SetAConst(dt);
    //ProductCoefficient coeff(dt, *r_alpha_dt);
    t.AddDomainIntegrator(new DiffusionIntegrator(*r_alpha_dt));
    t.Assemble();
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "aja..." << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    t.FormSystemMatrix(ess_tdof_list, T_op);

    if(T_prec) delete T_prec;
    T_prec = new HypreSmoother;
    T_solver.SetPreconditioner(*T_prec);
    T_solver.SetOperator(T_op);

    k->Mult(tmp_X, tmp_z);
    tmp_z.Neg();

    OperatorHandle A;
    Vector B;
    ParGridFunction tmp_dX_dt(&fespace);
    tmp_dX_dt = 0.0;

    t.FormLinearSystem(ess_tdof_list, tmp_dX_dt, tmp_z, A, dX_dt, B);

    T_solver.Mult(B, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &b,
                             int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    //if (T) delete T;
    //T = Add(1., M, scaled_dt, K);
    //T_solver.SetOperator(*T);
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
