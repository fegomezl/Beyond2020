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

        //Update information
        x->SetFromTrueDofs(*X);

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
             << progress << "\r";
        cout.flush();
    }
}

void Transport_Operator::SetParameters(const Vector &X){
    //Create the auxiliar grid functions
    aux.SetFromTrueDofs(X);

    //Associate the values of each auxiliar function
    for (int ii = 0; ii < aux.Size(); ii++){
        //Map option
        double H = aux(ii);
        aux(ii) = config.a_s*(0.5*(1 - tanh(5*JumpScale*(H+1))))
                + config.a_l*(0.5*(1 + tanh(5*JumpScale*(H-1))));

        //Conditional option
        /*
        if (aux(ii) >= 1)
            aux(ii) = config.a_l;
        else if (aux(ii) <= -1)
            aux(ii) = config.a_s;
        else
            aux(ii) = 0.;
        */
    }

    //Set the grid coefficients
    aux.ExchangeFaceNbrData();
    GridFunctionCoefficient coeff_A(&aux);
    coeff_rA.SetBCoef(coeff_A);
    
    config.sigma = -1.;
    config.kappa = pow((config.order+1),2);
    config.eta   = 0.;

    //Create RHS
    ConstantCoefficient bdr_s(H_s), bdr_l(H_l);
    if (b) delete b;
    if (B) delete B;
    b = new ParLinearForm(&fespace);
    b->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(bdr_l, coeff_rA_l, config.sigma, config.kappa), ess_bdr);
    //b->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(bdr_l, coeff_rA_l, config.sigma, config.kappa), ess_bdr_l);
    //b->AddBdrFaceIntegrator(new DGDirichletLFIntegrator(bdr_s, coeff_rA_s, config.sigma, config.kappa), ess_bdr_s);
    b->Assemble();
    B = b->ParallelAssemble();

    //Create transport matrix
    if (k) delete k;
    if (K) delete K;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rA_l));
    k->AddInteriorFaceIntegrator(new DGDiffusionIntegrator(coeff_rA_l, config.sigma, config.kappa));
    k->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA_l, config.sigma, config.kappa), ess_bdr);
    //k->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA_l, config.sigma, config.kappa), ess_bdr_l);
    //k->AddBdrFaceIntegrator(new DGDiffusionIntegrator(coeff_rA_s, config.sigma, config.kappa), ess_bdr_s);

    //Dont know how it works
    //k->AddInteriorFaceIntegrator(new DGDiffusionBR2Integrator(&fespace, config.eta));
    //k->AddBdrFaceIntegrator(new DGDiffusionBR2Integrator(&fespace, config.eta), ess_bdr);

    k->Assemble();
    k->Finalize();
    K = k->ParallelAssemble();    
}

void Transport_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //From  M(dX_dt) + K(X) = B
    //Solve M(dX_dt) + K(X) = B for dX_dt
    Z = 0.;
    dX_dt = 0.;
    
    K->Mult(-1., X, 1., Z);
    Z.Add(1., *B);

    M_solver.Mult(Z, dX_dt);
}

int Transport_Operator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + dt*K
    if (T) delete T;
    T = Add(1., *M, scaled_dt, *K);
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);

    //Set dt for RHS
    if (B_dt) delete B_dt;
    B_dt = new HypreParVector(&fespace);
    B_dt->Set(scaled_dt, *B);

    *j_status = 1;
    return 0;
}

int Transport_Operator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    //From  M(dX_dt) + K(X) = B
    //Solve M(X_new - X) + dt*K(X_new) = dt*B for X_new
    Z = 0.;
    X_new = X;

    M->Mult(X, Z);
    Z.Add(1., *B_dt);

    T_solver.Mult(Z, X_new);
    return 0;
}
