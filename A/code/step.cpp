#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt,  config.t_final - t);//****************************************    

    //Perform the time_step
    ode_solver->Step(X, t, dt);//******************************************
    oper->SetParameters(X);

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
    T = Add(1.0, M, dt, K);
    T_solver.SetOperator(*T);
    K.Mult(X, z);
    z.Neg();
    T_solver.Mult(z, dX_dt);
}
//***************************************************************************************************************
int Conduction_Operator::SUNImplicitSetup(const Vector &x,
                                         const Vector &fx, int jok, int *jcur,
                                         double gamma)
{
   // Setup the ODE Jacobian T = M + gamma K.
   if (T) delete T;
   T = Add(1., M, gamma, K);
   T_solver.SetOperator(*T);
   *jcur = 1;
   return (0);
}

int Conduction_Operator::SUNImplicitSolve(const Vector &b, Vector &x, double tol)
{
   // Solve the system A x = z => (M - gamma K) x = M b.
   M.Mult(b, z);
   T_solver.Mult(z, x);
   return (0);
}

