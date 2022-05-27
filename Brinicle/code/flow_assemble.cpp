#include "header.h"

//Boundary conditions
double boundary_vorticity_f(const Vector &x);

double boundary_stream_f(const Vector &x);
double boundary_stream_in_f(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1):
    config(config),
    fespace_H1(fespace_H1),
    block_offsets_H1(block_offsets_H1),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    ess_bdr_in(attributes),
    ess_bdr_closed_down(attributes), ess_bdr_closed_up(attributes),
    vorticity_boundary(&fespace_H1), stream_boundary(&fespace_H1), 
    velocity(&fespace_ND), rvelocity(&fespace_ND),
    stream(&fespace_H1), stream_gradient(&fespace_ND), 
    temperature(&fespace_H1), salinity(&fespace_H1), 
    density(&fespace_H1), density_dr(&fespace_H1), 
    impermeability(&fespace_H1), 
    B(block_offsets_H1),
    B0(NULL), B1(NULL),
    A00(NULL), A01(NULL), A10(NULL), A11(NULL),
    coeff_r(r_f), coeff_r_inv(r_inv_f), 
    coeff_r_inv_hat(dim, r_inv_hat_f),
    coeff_rot(dim, rot_f), 
    coeff_vorticity(boundary_vorticity_f), 
    coeff_stream(boundary_stream_f), 
    coeff_stream_in(boundary_stream_in_f),
    coeff_stream_closed_down(0.), coeff_stream_closed_up(InflowFlux),
    gradient(&fespace_H1, &fespace_ND)
{ 
    /****
     * Define essential boundary conditions
     * 
     *           4              1
     *         /---|---------------------------\
     *         |                               |
     *         |                               |
     *         |                               | 
     *         |                               |
     *       2 |                               | 3
     *         |                               |
     *         |                               |
     *         |                               | 
     *         |                               |
     *         \-------------------------------/
     *                          0
     *
     ****/
    
    //Vorticity boundary conditions
    ess_bdr_0 = 0;
    ess_bdr_0[0] = 0; ess_bdr_0[1] = 0;
    ess_bdr_0[2] = 1; ess_bdr_0[3] = 0;
    ess_bdr_0[4] = 0; 
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);
  
    //Stream boundary conditions
    ess_bdr_1 = 0;
    ess_bdr_1[0] = 0; ess_bdr_1[1] = 1;
    ess_bdr_1[2] = 1; ess_bdr_1[3] = 1;
    ess_bdr_1[4] = 1; 
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    ess_bdr_in = 0;
    ess_bdr_in[4] = 1;
     
    ess_bdr_closed_down = 0;
    ess_bdr_closed_down[0] = 1;
    ess_bdr_closed_down[2] = 1;
    
    ess_bdr_closed_up = 0;
    ess_bdr_closed_up[1] = 1;
    ess_bdr_closed_up[3] = 1;
    ess_bdr_closed_up[5] = 1;

    //Apply boundary conditions
    vorticity_boundary.ProjectCoefficient(coeff_vorticity);
    vorticity_boundary.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);

    stream_boundary.ProjectCoefficient(coeff_stream);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_in);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_closed_down);
    stream_boundary.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_closed_up);

    //Define the constant RHS
    ParLinearForm b0(&fespace_H1);
    ConstantCoefficient Zero(0.);
    b0.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    b0.Assemble();

    //Define constant bilinear forms of the system
    ParBilinearForm a00(&fespace_H1);
    a00.AddDomainIntegrator(new MassIntegrator);
    a00.Assemble();
    a00.EliminateEssentialBC(ess_bdr_0, vorticity_boundary, b0, Operator::DIAG_ONE);
    a00.Finalize();
    A00 = a00.ParallelAssemble();

    ParMixedBilinearForm a01(&fespace_H1, &fespace_H1);
    a01.AddDomainIntegrator(new MixedGradGradIntegrator);
    a01.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a01.Assemble();
    a01.EliminateTrialDofs(ess_bdr_1, stream_boundary, b0);
    a01.EliminateTestDofs(ess_bdr_0);    
    a01.Finalize();
    A01 = a01.ParallelAssemble();

    B0 = b0.ParallelAssemble();

    //Create gradient interpolator
    gradient.AddDomainIntegrator(new GradientInterpolator);
    gradient.Assemble();
    gradient.Finalize();
}

//Boundary condition for vorticity
double boundary_vorticity_f(const Vector &x){
    return 0.;
}

//Boundary condition for stream
double boundary_stream_f(const Vector &x){
    double x_rel = x(0)/RInflow;
    double y_rel = x(1)/Z;
    if (x(0) < RInflow)
        return InflowFlux*y_rel*pow(x_rel, 2)*(2-pow(x_rel, 2));
    return InflowFlux*y_rel;
}

//Boundary condition for stream at inlet
double boundary_stream_in_f(const Vector &x){
    double x_rel = x(0)/RInflow;
    return InflowFlux*pow(x_rel, 2)*(2-pow(x_rel, 2));
}
