#include "header.h"

//Boundary values
double bdr_vorticity(const Vector &x);
double bdr_stream(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1):
    config(config),
    fespace_H1(fespace_H1), fespace_L2(fespace_L2), fespace_ND(fespace_ND),
    block_offsets_H1(block_offsets_H1),
    vorticity(&fespace_H1), stream(&fespace_H1),
    Vorticity(&fespace_H1), Stream(&fespace_H1),
    A00(NULL),   A01(NULL),   A10(NULL),   A11(NULL),
    A00_e(NULL), A01_e(NULL), A10_e(NULL), A11_e(NULL),
    B0(&fespace_H1), B1(&fespace_H1),
    coeff_r(r_f), coeff_r_inv(r_inv_f), 
    coeff_r_inv_hat(dim, r_inv_hat_f), 
    coeff_rot(dim, rot_f), 
    gradient(&fespace_H1, &fespace_ND)
{  
    //Set boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //           2|            |3
    //            |            |
    //            \------------/
    //                  0
    //

    Array<int> ess_bdr_0(attributes);   
    ess_bdr_0 = 0;
    ess_bdr_0[0] = 1; ess_bdr_0[1] = 1;
    ess_bdr_0[2] = 1; ess_bdr_0[3] = 1;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);
  
    Array<int> ess_bdr_1(attributes);   
    ess_bdr_1 = 0;
    ess_bdr_1[0] = 1; ess_bdr_1[1] = 1;
    ess_bdr_1[2] = 1; ess_bdr_1[3] = 1;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    //Setup boundary coefficients
    FunctionCoefficient coeff_bdr_vorticity(bdr_vorticity);
    FunctionCoefficient coeff_bdr_stream(bdr_stream);

    //Apply boundary conditions
    vorticity.ProjectCoefficient(coeff_bdr_vorticity);
    vorticity.ParallelProject(Vorticity);
 
    stream.ProjectCoefficient(coeff_bdr_stream);
    stream.ParallelProject(Stream);

    //Define constant bilinear forms of the system
    ParBilinearForm a00(&fespace_H1);
    a00.AddDomainIntegrator(new MassIntegrator);
    a00.Assemble();
    a00.Finalize();
    A00 = a00.ParallelAssemble();
    A00_e = A00->EliminateRowsCols(ess_tdof_0);

    ParMixedBilinearForm a01(&fespace_H1, &fespace_H1);
    a01.AddDomainIntegrator(new MixedGradGradIntegrator);
    a01.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a01.Assemble();
    a01.Finalize();
    A01 = a01.ParallelAssemble();
    A01->EliminateRows(ess_tdof_0);
    A01_e = A01->EliminateCols(ess_tdof_1);

    A10 = A01->Transpose();
    A10->EliminateRows(ess_tdof_1);
    A10_e = A10->EliminateCols(ess_tdof_0);

    //Define the constant RHS
    ParLinearForm b0(&fespace_H1);
    b0.Assemble();
    b0.ParallelAssemble(B0);
    A01_e->Mult(Stream, B0, -1., 1.);
    EliminateBC(*A00, *A00_e, ess_tdof_0, Vorticity, B0);

    //Create gradient interpolator
    gradient.AddDomainIntegrator(new GradientInterpolator);
    gradient.Assemble();
    gradient.Finalize();
}

double bdr_vorticity(const Vector &x){
    return 0.;
}

double bdr_stream(const Vector &x){
    return 0.;
}
