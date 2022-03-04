#include "header.h"

//Boundary values
double boundary_vorticity_f(const Vector &x);

double boundary_stream_f(const Vector &x);
double boundary_stream_in_f(const Vector &x);
double boundary_stream_out_f(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X):
    config(config),
    fespace_H1(fespace_H1),
    block_offsets_H1(block_offsets_H1),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    Vorticity(&fespace_H1), Stream(&fespace_H1),
    A00(NULL),  A00_e(NULL), 
    A01(NULL),  A01_e(NULL), 
    A10(NULL),  A10_e(NULL), 
    A11(NULL),  A11_e(NULL),
    B0(NULL), B1(NULL),
    B(block_offsets_H1),
    velocity(&fespace_ND), rvelocity(&fespace_ND),
    temperature(&fespace_H1), temperature_dr(&fespace_H1), 
    salinity(&fespace_H1), salinity_dr(&fespace_H1), 
    phase(&fespace_H1), impermeability(&fespace_H1), 
    stream(&fespace_H1), stream_gradient(&fespace_ND), 
    coeff_r(r_f), coeff_r_inv(r_inv_f), coeff_r_inv_hat(dim, r_inv_hat_f),
    gradient(&fespace_H1, &fespace_ND),
    coeff_rot(dim, rot_f)
{ 
    //Define essential boundary conditions
    //   
    //              4               1
    //            /---|---------------------------\
    //            |                               |
    //            |                               |
    //            |                               | 3
    //            |                               |
    //          2 |                               |
    //            |                               -
    //            |                               |
    //            |                               | 5
    //            |                               |
    //            \-------------------------------/
    //                            0
    //

    ess_bdr_0 = 0;
    ess_bdr_0[0] = 0; ess_bdr_0[1] = 0;
    ess_bdr_0[2] = 1; ess_bdr_0[3] = 0;
    ess_bdr_0[4] = 0; ess_bdr_0[5] = 0;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);
  
    ess_bdr_1 = 0;
    ess_bdr_1[0] = 1; ess_bdr_1[1] = 1;
    ess_bdr_1[2] = 1; ess_bdr_1[3] = 1;
    ess_bdr_1[4] = 1; ess_bdr_1[5] = 1;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    Array<int> ess_bdr_1_in(attributes);
    ess_bdr_1_in = 0;
    ess_bdr_1_in[4] = ess_bdr_1[4];
    
    Array<int> ess_bdr_1_out(attributes);
    ess_bdr_1_out = 0;
    ess_bdr_1_out[5] = ess_bdr_1[5];
    
    Array<int> ess_bdr_1_closed_down(attributes);
    ess_bdr_1_closed_down = 0;
    ess_bdr_1_closed_down[0] = ess_bdr_1[0];
    ess_bdr_1_closed_down[2] = ess_bdr_1[2];
    
    Array<int> ess_bdr_1_closed_up(attributes);
    ess_bdr_1_closed_up = 0;
    ess_bdr_1_closed_up[1] = ess_bdr_1[1];
    ess_bdr_1_closed_up[3] = ess_bdr_1[3];

    //Apply boundary conditions
    FunctionCoefficient coeff_vorticity(boundary_vorticity_f);
    ParGridFunction vorticity(&fespace_H1);
    vorticity.ProjectCoefficient(coeff_vorticity);
    vorticity.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);
    vorticity.ParallelProject(Vorticity);
 
    FunctionCoefficient coeff_stream(boundary_stream_f);
    FunctionCoefficient coeff_stream_in(boundary_stream_in_f);
    FunctionCoefficient coeff_stream_out(boundary_stream_out_f);
    ConstantCoefficient coeff_stream_closed_down(0.);
    ConstantCoefficient coeff_stream_closed_up(InflowFlux);
    ParGridFunction stream(&fespace_H1);
    stream.ProjectCoefficient(coeff_stream);
    stream.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_1_in);
    stream.ProjectBdrCoefficient(coeff_stream_out, ess_bdr_1_out);
    stream.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_1_closed_down);
    stream.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_1_closed_up);
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
    ConstantCoefficient Zero(0.);
    b0.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    b0.Assemble();
    B0 = b0.ParallelAssemble();
    A10_e->Mult(Stream, *B0, -1., 1.);
    EliminateBC(*A00, *A00_e, ess_tdof_0, Vorticity, *B0);

    //Create gradient interpolator
    gradient.AddDomainIntegrator(new GradientInterpolator);
    gradient.Assemble();
    gradient.Finalize();
}

//Boundary values
double boundary_vorticity_f(const Vector &x){
    return 0.;
}

double boundary_stream_f(const Vector &x){
    double x_rel = x(0)/RIn;
    double y_rel = x(1)/ZOut;
    double in = 1., out = 1.;
    if (x(0) < RIn)
        in = pow(x_rel, 2)*(2-pow(x_rel, 2));
    if (x(1) < ZOut)
        out = pow(y_rel, 2)*(3-2*y_rel);
    return InflowFlux*in*out;
}

double boundary_stream_in_f(const Vector &x){
    double x_rel = x(0)/RIn;
    return InflowFlux*pow(x_rel, 2)*(2-pow(x_rel, 2));
}

double boundary_stream_out_f(const Vector &x){
    double y_rel = x(1)/ZOut;
    return InflowFlux*pow(y_rel, 2)*(3-2*y_rel);
}
