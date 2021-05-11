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
    double c_s = 1.878, c_l = 4.219;
    double k_s = 133.2, k_l = 33.6;
    double L = 302.3;

    double DeltaT = 0.0001;

    Vector Aux(X);
    if (T_f != 0)
        Aux -= T_f;

    //Create the auxiliar grid functions
    ParGridFunction aux(&fespace);   aux.SetFromTrueDofs(Aux);
    ParGridFunction aux_C(&fespace); aux.SetFromTrueDofs(Aux);
    ParGridFunction aux_K(&fespace); aux.SetFromTrueDofs(Aux);

    //Associated the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        if (aux(ii) > 0){
            aux_C(ii) = c_l;
            aux_K(ii) = k_l;
        } else {
            aux_C(ii) = c_s;
            aux_K(ii) = k_s;
        }

        aux(ii) = (L/(DeltaT*sqrt(2*M_PI)))*exp(-pow(aux(ii)/DeltaT, 2)/2);
    }

    //Create the associated coefficients
    GridFunctionCoefficient coeff_C(&aux_C); ProductCoefficient coeff_rC(r, coeff_C);
    GridFunctionCoefficient coeff_K(&aux_K); ProductCoefficient coeff_rK(r, coeff_K);
    GridFunctionCoefficient coeff_L(&aux);   ProductCoefficient coeff_rL(r, coeff_L);

    //Create corresponding bilinear forms
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rC));
    m->AddDomainIntegrator(new MassIntegrator(coeff_rL));
    m->Assemble(0);
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->Assemble(0);
    k->FormSystemMatrix(ess_tdof_list, K);
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
