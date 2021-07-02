#include "header.h"

void Artic_sea::time_step(){
    //Update iteration parameters
    last = (t >= config.t_final - 1e-8*config.dt_init);
    dt = min(dt, config.t_final - t);

    //Perform the time_step
    oper_T->SetParameters(*X_T);
    ode_solver->Step(*X_T, t, dt);

    //flow_oper->Solve(config, *x_T, X_Psi, x_psi);

    //Update visualization steps
    vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);

    if (last || vis_steps <= vis_iteration){
        //Update parameters
        vis_iteration = 0;
        vis_impressions++;

        //Graph
        x_T->SetFromTrueDofs(*X_T);
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
    GridFunctionCoefficient coeff_C(&aux_C);
    GridFunctionCoefficient coeff_K(&aux_K);
    GridFunctionCoefficient coeff_L(&aux);

    SumCoefficient coeff_CL(coeff_C, coeff_L);

    coeff_rK.SetBCoef(coeff_K); dt_coeff_rK.SetBCoef(coeff_rK); 
    coeff_rCL.SetBCoef(coeff_CL);
    coeff_rCLV.SetACoef(coeff_CL);
    dt_coeff_rCLV.SetBCoef(coeff_rCLV);

    //Create corresponding bilinear forms
    delete m;
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    m->Assemble();
    m->FormSystemMatrix(ess_tdof_list, M);
    M_solver.SetOperator(M);

    delete k;
    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(coeff_rK));
    k->AddDomainIntegrator(new ConvectionIntegrator(coeff_rCLV));
    k->Assemble();
    k->Finalize();
}

void Conduction_Operator::UpdateVelocity(const HypreParVector &Psi){
    psi.SetFromTrueDofs(Psi);
    gradpsi.SetGridFunction(&psi);
    rV.SetBCoef(gradpsi);
    coeff_rCLV.SetBCoef(rV);
}

void Flow_Operator::Solve(Config config, const ParGridFunction &x_T, HypreParVector *X_Psi, ParGridFunction *x_psi){
    //Create preconditioner objects
    HypreParVector *Dd = new HypreParVector(MPI_COMM_WORLD, D->GetGlobalNumRows(),
                                            D->GetRowStarts());
    //HypreParVector *Md = new HypreParVector(MPI_COMM_WORLD, M->GetGlobalNumRows(),
    //                                        M->GetRowStarts());
    HypreParMatrix *Dd_inv_Ct = NULL;
    HypreParMatrix *C_Dd_inv_Ct = NULL;
    HypreParMatrix *S = NULL;
    HypreParMatrix *Mplus = NULL;

    //M->GetDiag(*Md);
    D->GetDiag(*Dd);
    //Md_inv_Ct = C->Transpose();
    Dd_inv_Ct = Ct;
    Dd_inv_Ct->InvScaleRows(*Dd);
    C_Dd_inv_Ct = ParMult(C, Dd_inv_Ct);              //S = +/-(M - C diag(D)^-1 C^t)
    Mplus = M;
    (*Mplus) *= -1.;
    S = ParAdd(Mplus, C_Dd_inv_Ct);

    //Create solvers for the preconditioner
    Solver *D_inv = new HypreDiagScale(*D);
    D_inv->iterative_mode = false;

    Solver *S_inv = new HypreBoomerAMG(*S);
    S_inv->iterative_mode = false;

    //Create the complete preconditioner
    //
    //   P = [ diag(D)  0]
    //       [   0      S]
    BlockDiagonalPreconditioner *Prec = new BlockDiagonalPreconditioner(block_true_offsets);
    Prec->SetDiagonalBlock(0, D_inv);
    Prec->SetDiagonalBlock(1, S_inv);

    //Solve the linear system Ax=B
    MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetPreconditioner(*Prec);
    solver.SetOperator(*A);
    solver.SetPrintLevel(0);
    solver.SetAbsTol(1e-10);
    solver.SetRelTol(1e-6);
    solver.SetMaxIter(1000);
    Y = 0.;
    solver.Mult(B, Y);

    //Recover the solution on each proccesor
    psi->MakeRef(&fespace, y.GetBlock(0), 0);
    psi->Distribute(&(Y.GetBlock(0)));

    w->MakeRef(&fespace, y.GetBlock(1), 0);
    w->Distribute(&(Y.GetBlock(1)));

    GradientGridFunctionCoefficient psi_grad(psi);

    X_Psi = psi->GetTrueDofs();
    x_psi->SetFromTrueDofs(*X_Psi);

    //Delete used memory
    delete Dd, Mplus, Dd_inv_Ct, C_Dd_inv_Ct, S;
    delete D_inv, S_inv;
    delete Prec;
}
