#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Create the block offsets
    Array<int> block_true_offsets(3);
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace.TrueVSize();
    block_true_offsets[2] = fespace.TrueVSize();
    block_true_offsets.PartialSum();

    //Initialize the corresponding vectors
    HypreParVector X_aux(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        X_aux(ii) = X(ii);

    //Solve M(dX_dt) = -K(X) for dX_dt
    ParGridFunction x(&fespace);
    ParLinearForm z(&fespace);
    x.SetFromTrueDofs(X_aux);
    k->Mult(x, z);
    z.Neg();

    OperatorHandle A;
    HypreParVector Z(&fespace);
    ParGridFunction dx_dt(&fespace);
    HypreParVector dX_dt_aux(&fespace);
    dx_dt = 0.;

    m->FormLinearSystem(ess_tdof_list, dx_dt, z, A, dX_dt_aux, Z);
    M_solver.Mult(Z, dX_dt_aux);

    dX_dt = 0.;
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        dX_dt(ii) = dX_dt_aux(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        dX_dt(ii) = dX_dt_aux(ii - block_true_offsets[1]);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    if (t) delete t;
    t = new ParBilinearForm(&fespace);
    dt_coeff_rK.SetAConst(scaled_dt); dt_coeff_rCLV.SetAConst(scaled_dt);
    t->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    t->AddDomainIntegrator(new ConvectionIntegrator(dt_coeff_rCLV));
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rK));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T);
    T_solver.SetOperator(T);

    *j_status = 1;
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &B, Vector &X, double tol){
    //Solve the system Ax = z -> (M + gamma*K)x = Mb
    
    //Create the block offsets
    Array<int> block_true_offsets(3);
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace.TrueVSize();
    block_true_offsets[2] = fespace.TrueVSize();
    block_true_offsets.PartialSum();

    //Initialize the corresponding vectors
    HypreParVector B_aux(&fespace);
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        B_aux(ii) = B(ii);
    HypreParVector X_aux(&fespace);

    HypreParVector Z(&fespace);
    M.Mult(B_aux, Z);
    T_solver.Mult(Z, X_aux);

    X = 0.;
    for (int ii = block_true_offsets[0]; ii < block_true_offsets[1]; ii++)
        X(ii) = X_aux(ii);
    for (int ii = block_true_offsets[1]; ii < block_true_offsets[2]; ii++)
        X(ii) = X_aux(ii - block_true_offsets[1]);
    return 0;
}
