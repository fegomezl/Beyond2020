#include "header.h"

void Artic_sea::time_step(){

  //adaptative mesh***************************************************************************
  refiner.Reset();
  derefiner.Reset();
  refiner.Apply(pmesh);
  // 21. Quit the AMR loop if the termination criterion has been met
  if (refiner.Stop())
    {
      a.Update(); // Free the assembled data
      break;
    }
  UpdateAndRebalance(pmesh, fespace, x, a, b);
  if (derefiner.Apply(pmesh))
    {
      // 24. Update the space and the solution, rebalance the mesh.
      UpdateAndRebalance(pmesh, fespace, x, a, b);
    }
  //*******************************************************************************************
  
  //Check for last iteration
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
	
        //Calculate convergence
        initial_f.SetTime(t);
        x->SetFromTrueDofs(*X);
        actual_error = x->ComputeL2Error(initial_f);
        actual_error /= (Zmax - Zmin)*(Rmax - Rmin);
        total_error += actual_error;

        //Graph
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
             << progress << setw(9)
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
    ProductCoefficient coeff(alpha, r);
    
    delete k;
    k = new ParBilinearForm(&fespace);
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
//*******************************************************************************
void UpdateAndRebalance(ParMesh &pmesh, ParFiniteElementSpace &fespace,
                        ParGridFunction &x, ParBilinearForm &a,
                        ParLinearForm &b)
{
   // Update the space: recalculate the number of DOFs and construct a matrix
   // that will adjust any GridFunctions to the new mesh state.
   fespace.Update();

   // Interpolate the solution on the new mesh by applying the transformation
   // matrix computed in the finite element space. Multiple GridFunctions could
   // be updated here.
   x.Update();

   if (pmesh.Nonconforming())
   {
      // Load balance the mesh.
      pmesh.Rebalance();

      // Update the space again, this time a GridFunction redistribution matrix
      // is created. Apply it to the solution.
      fespace.Update();
      x.Update();
   }

   // Inform the linear and bilinear forms that the space has changed.
   a.Update();
   b.Update();

   // Free any transformation matrices to save memory.
   fespace.UpdatesFinished();
}
