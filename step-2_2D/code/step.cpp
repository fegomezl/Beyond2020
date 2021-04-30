#include "header.h"

void Artic_sea::time_step(){
  //Check for last iteration
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    oper->SetParameters(*X);
    ode_solver->Step(*X, t, dt);

    double actual_error;
    if (last || (iteration % config.vis_steps) == 0){
      //Calculate convergence
      initial_f.SetTime(t);
      x->SetFromTrueDofs(*X);
      actual_error = x->ComputeL2Error(initial_f);

      //Graph
      paraview_out->SetCycle(iteration);
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
             << t  << setw(12)
             << progress << setw(12)
             << actual_error << "\r";
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
    ProductCoefficient coeff1(alpha, r);
    ScalarVectorProductCoefficient coeff2(alpha, r_hat);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff1));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff2));
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
