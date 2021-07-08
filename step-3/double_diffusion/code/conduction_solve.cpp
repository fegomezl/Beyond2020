#include "header.h"

void Conduction_Operator::Mult(const Vector &X, Vector &dX_dt) const{
    //Solve M(dX_dt) = -K(X) for dX_dt
    BlockVector Z(block_true_offsets);
    K->Mult(X, Z);
    Z.Neg();
    M_solver.Mult(Z, dX_dt);
}

int Conduction_Operator::SUNImplicitSetup(const Vector &X, const Vector &B, int j_update, int *j_status, double scaled_dt){
    //Setup the ODE Jacobian T = M + gamma*K
    /*if (t) delete t;
    t = new ParBilinearForm(&fespace);
    dt_coeff_rK.SetAConst(scaled_dt); dt_coeff_rCLV.SetAConst(scaled_dt);
    t->AddDomainIntegrator(new MassIntegrator(coeff_rCL));
    t->AddDomainIntegrator(new ConvectionIntegrator(dt_coeff_rCLV));
    t->AddDomainIntegrator(new DiffusionIntegrator(dt_coeff_rK));
    t->Assemble();
    t->FormSystemMatrix(ess_tdof_list, T);
    T_solver.SetOperator(T);

    *j_status = 1;*/
    return 0;
}

int Conduction_Operator::SUNImplicitSolve(const Vector &B, Vector &X, double tol){
    //Solve the system Ax = z -> (M - gamma*K)x = Mb
    /*HypreParVector Z(&fespace);
    M.Mult(B,Z);
    T_solver.Mult(Z,X);*/
    return 0;
}
