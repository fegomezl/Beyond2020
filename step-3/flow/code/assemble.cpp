#include "header.h"

//Right hand side of the equation
void f_rhs(const Vector &x, Vector &f){                 
    double r = x(0);
    double z = x(1);
    f(0) = z*pow(r, 2);
    f(1) = z*r*(height - z);
}

//Constant to the brinkman term
double porous_constant(const Vector &x){
    return x(1);
}

//Exact solution for the pressure
double p_exact(const Vector &x){               
    return (height/6)*pow(x(0), 3);
}

//Exact solution for the velocity
void v_exact(const Vector &x, Vector &f){
    double r = x(0);
    double z = x(1);
    f(0) = r*(z - height/2);
    f(1) = z*(height - z);
}

double r_f(const Vector &x){
    return x(0);
}

void r_hat(const Vector &x, Vector &f){
    f(0) = 1;
    f(1) = 0;
}

void Artic_sea::assemble_system(){
    //Set the boundary values
    Array<int> ess_tdof_list;

    //Define local coefficients
    VectorFunctionCoefficient f_coeff(dim, f_rhs);
    ConstantCoefficient g_coeff(0.);
    FunctionCoefficient eta_coeff(porous_constant);

    //Calculate cylindrical coefficients
    FunctionCoefficient r_coeff(r_f);
    ProductCoefficient neg_r_coeff(-1., r_coeff);
    VectorFunctionCoefficient r_hat_coeff(dim, r_hat);
    ScalarVectorProductCoefficient f_r(r_coeff, f_coeff);
    ProductCoefficient eta_r(r_coeff, eta_coeff);
    ProductCoefficient p_exact_r(neg_r_coeff, p_exact_coeff);

    //Define the RHS
    f = new ParLinearForm;
    f->Update(fespace_rt, b.GetBlock(0), 0);
    f->AddDomainIntegrator(new VectorFEDomainLFIntegrator(f_r));
    f->AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator(p_exact_r));
    f->Assemble();
    f->SyncAliasMemory(b);
    f->ParallelAssemble(B.GetBlock(0));
    B.GetBlock(0).SyncAliasMemory(B);

    g = new ParLinearForm;
    g->Update(fespace_l2, b.GetBlock(1), 0);
    g->AddDomainIntegrator(new DomainLFIntegrator(g_coeff));
    g->Assemble();
    g->SyncAliasMemory(b);
    g->ParallelAssemble(B.GetBlock(1));
    B.GetBlock(1).SyncAliasMemory(B);

    //Define bilinear forms of the system
    m = new ParBilinearForm(fespace_rt);
    m->AddDomainIntegrator(new VectorFEMassIntegrator(eta_r));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();

    c = new ParMixedBilinearForm(fespace_rt, fespace_l2);
    c->AddDomainIntegrator(new VectorFEDivergenceIntegrator(r_coeff));
    c->AddDomainIntegrator(new MixedDotProductIntegrator(r_hat_coeff));
    c->Assemble();
    c->Finalize();
    C = c->ParallelAssemble();
    (*C) *= -1;

    Ct = new TransposeOperator(C);

    //Create the complete bilinear operator:
    //
    //   A = [ M  C^t]
    //       [ C   0 ] 
    A = new BlockOperator(block_true_offsets);
    A->SetBlock(0, 0, M);
    A->SetBlock(0, 1, Ct);
    A->SetBlock(1, 0, C);
}
