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
    int n = 2;
    double DeltaT = 10;
    double inv_alpha_l = 0.126, inv_alpha_s = 0.014;
    double L_l = 0.3589, L_s = 0.0906;

    ParGridFunction aux(&fespace);
    aux.SetFromTrueDofs(X);

    //Create the K coefficient
    for (int ii = 0; ii < aux.Size(); ii++){
        double ref = aux(ii) - T_f;
        if (ref > 0){
            aux(ii) = inv_alpha_l;
            if (abs(ref) <= n*DeltaT)
                aux(ii) += L_l*exp(-pow(ref/DeltaT, 2)/2);
        } else {
            aux(ii) = inv_alpha_s;
            if (abs(ref) <= n*DeltaT)
                aux(ii) += L_s*exp(-pow(ref/DeltaT, 2)/2);
        }
    }
    GridFunctionCoefficient alpha(&aux);
    ProductCoefficient coeff(alpha, r);

    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff));
    m->Assemble(0);
    m->FormSystemMatrix(ess_tdof_list, M);
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
