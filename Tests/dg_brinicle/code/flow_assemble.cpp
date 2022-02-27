#include "header.h"

//Boundary values
double boundary_vorticity_f(const Vector &x);
double boundary_stream_f(const Vector &x);
double boundary_stream_in_f(const Vector &x);
double boundary_stream_out_f(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1):
    config(config),
    fespace_H1(fespace_H1), fespace_L2(fespace_L2), fespace_ND(fespace_ND),
    block_offsets_H1(block_offsets_H1),
    vorticity(&fespace_H1), stream(&fespace_H1),
    Vorticity(&fespace_H1), Stream(&fespace_H1),
    A00(NULL),  A00_e(NULL), 
    A01(NULL),  A01_e(NULL), 
    A10(NULL),  A10_e(NULL), 
    A11(NULL),  A11_e(NULL),
    B0(&fespace_H1), B1(&fespace_H1),
    coeff_r(r_f), coeff_r_inv(r_inv_f), 
    coeff_r_inv_hat(dim, r_inv_hat_f), 
    coeff_rot(dim, rot_f), 
    gradient(&fespace_H1, &fespace_ND)
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
    Array<int> ess_bdr_0(attributes);   
    ess_bdr_0 = 0;
    ess_bdr_0[0] = 1; ess_bdr_0[1] = 1;
    ess_bdr_0[2] = 1; ess_bdr_0[3] = 1;
    ess_bdr_0[4] = 1; ess_bdr_0[5] = 1;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);
  
    Array<int> ess_bdr_1(attributes);   
    ess_bdr_1 = 0;
    ess_bdr_1[0] = 1; ess_bdr_1[1] = 1;
    ess_bdr_1[2] = 1; ess_bdr_1[3] = 1;
    ess_bdr_1[3] = 1; ess_bdr_1[5] = 1;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    Array<int> bdr_stream_in(attributes);   
    bdr_stream_in = 0;
    bdr_stream_in[4] = ess_bdr_1[4];
    
    Array<int> bdr_stream_out(attributes);   
    bdr_stream_out = 0;
    bdr_stream_out[5] = ess_bdr_1[5];
    
    Array<int> bdr_stream_closed_down(attributes);   
    bdr_stream_closed_down = 0;
    bdr_stream_closed_down[0] = ess_bdr_1[0];
    bdr_stream_closed_down[2] = ess_bdr_1[2];
    
    Array<int> bdr_stream_closed_up(attributes);   
    bdr_stream_closed_up = 0;
    bdr_stream_closed_up[1] = ess_bdr_1[1];
    bdr_stream_closed_up[3] = ess_bdr_1[3];

    //Setup boundary coefficients
    FunctionCoefficient vorticity_coeff(boundary_vorticity_f);
    FunctionCoefficient stream_coeff(boundary_stream_f);
    FunctionCoefficient stream_in_coeff(boundary_stream_in_f);
    FunctionCoefficient stream_out_coeff(boundary_stream_out_f);
    ConstantCoefficient stream_closed_up_coeff(InflowFlux);
    ConstantCoefficient stream_closed_down_coeff(0.);

    //Apply boundary conditions
    vorticity.ProjectCoefficient(vorticity_coeff);
    vorticity.ParallelProject(Vorticity);
 
    stream.ProjectCoefficient(stream_coeff);
    stream.ProjectBdrCoefficient(stream_in_coeff, bdr_stream_in);
    stream.ProjectBdrCoefficient(stream_out_coeff, bdr_stream_out);
    stream.ProjectBdrCoefficient(stream_closed_up_coeff, bdr_stream_closed_up);
    stream.ProjectBdrCoefficient(stream_closed_down_coeff, bdr_stream_closed_down);
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

double boundary_vorticity_f(const Vector &x){
    return 0.;
}

double boundary_stream_f(const Vector &x){
    double Q = 0.25*InflowVelocity*pow(R_in, -2);
    double x_rel = x(0)/R_in;
    double y_rel = x(1)/Z_out;
    double in = 1., out = 1.;
    if (x(0) < R_in)
        in = pow(x_rel, 2)*(2-pow(x_rel, 2));
    if (x(1) < Z_out)
        out = pow(y_rel, 2)*(3-2*y_rel);
    return Q*in*out;
}

double boundary_stream_in_f(const Vector &x){
    double Q = 0.25*InflowVelocity*pow(R_in, -2);
    double x_rel = x(0)/R_in;
    return Q*pow(x_rel, 2)*(2-pow(x_rel, 2));
}

double boundary_stream_out_f(const Vector &x){
    double Q = 0.25*InflowVelocity*pow(R_in, -2);
    double y_rel = x(1)/Z_out;
    return Q*pow(y_rel, 2)*(3-2*y_rel);
}
