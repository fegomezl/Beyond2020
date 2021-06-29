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
    //Create the auxiliar grid functions
    Vector Aux(X);
    if (config.T_f != 0)
        Aux -= config.T_f;
    aux.SetFromTrueDofs(Aux);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > 0){
            aux_C(ii) = config.c_l;
            aux_K(ii) = config.k_l;
        } else {
            aux_C(ii) = config.c_s;
            aux_K(ii) = config.k_s;
        }

        aux(ii) = config.L*config.invDeltaT*exp(-M_PI*pow(config.invDeltaT*aux(ii), 2));
    }

    //Set the associated coefficients
    coeff_C.SetGridFunction(&aux_C);
    coeff_K.SetGridFunction(&aux_K);
    coeff_L.SetGridFunction(&aux);

    coeff_rC.SetBCoef(coeff_C);
    coeff_rK.SetBCoef(coeff_K); 
    coeff_rL.SetBCoef(coeff_L);

    //Create corresponding bilinear forms
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m->AddDomainIntegrator(new MassIntegrator(coeff_rL));
    m->Assemble();
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->Assemble();
    k->Finalize();
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
    dt_coeff_rK.SetAConst(dt); dt_coeff_rK.SetBCoef(coeff_rK);
    t->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    t->AddDomainIntegrator(new MassIntegrator(coeff_rL));
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rK));
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
    dt_coeff_rK.SetAConst(scaled_dt); dt_coeff_rK.SetBCoef(coeff_rK);
    t->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    t->AddDomainIntegrator(new MassIntegrator(coeff_rL));
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rK));
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
