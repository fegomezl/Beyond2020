#include "header.h"

//Boundary conditions
double boundary_stream_f(const Vector &x);
double boundary_stream_in_f(const Vector &x);

Flow_Operator::Flow_Operator(Config config, ParFiniteElementSpace &fespace_H1, int dim, int attributes, Array<int> block_offsets_H1):
    config(config),
    fespace_H1(fespace_H1),
    block_offsets_H1(block_offsets_H1),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    ess_bdr_in(attributes),
    ess_bdr_closed_down(attributes), ess_bdr_closed_up(attributes),
    vorticity(&fespace_H1), stream(&fespace_H1), 
    temperature(&fespace_H1), salinity(&fespace_H1), 
    density(&fespace_H1), density_dr(&fespace_H1), 
    impermeability(&fespace_H1), 
    B(block_offsets_H1),
    B0(&fespace_H1), B1(&fespace_H1),
    A00(NULL), A01(NULL), A10(NULL), A11(NULL),
    coeff_r(r_f), coeff_r_inv_hat(dim, r_inv_hat_f),
    coeff_vorticity(0.), coeff_stream(boundary_stream_f),
    coeff_stream_in(boundary_stream_in_f),
    coeff_stream_closed_down(0.),
    coeff_stream_closed_up(1.)
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
    vorticity.ProjectCoefficient(coeff_vorticity);
    vorticity.ProjectBdrCoefficient(coeff_vorticity, ess_bdr_0);

    stream.ProjectCoefficient(coeff_stream);
    stream.ProjectBdrCoefficient(coeff_stream_in, ess_bdr_in);
    stream.ProjectBdrCoefficient(coeff_stream_closed_down, ess_bdr_closed_down);
    stream.ProjectBdrCoefficient(coeff_stream_closed_up, ess_bdr_closed_up);

    //Define the constant RHS
    ParLinearForm b0(&fespace_H1);
    ConstantCoefficient Zero(0.);
    b0.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    b0.Assemble();

    //Define constant bilinear forms of the system
    ParBilinearForm a00(&fespace_H1);
    a00.AddDomainIntegrator(new MassIntegrator);
    a00.Assemble();
    a00.EliminateEssentialBC(ess_bdr_0, vorticity, b0, Operator::DIAG_ONE);
    a00.Finalize();
    A00 = a00.ParallelAssemble();

    ParMixedBilinearForm a01(&fespace_H1, &fespace_H1);
    a01.AddDomainIntegrator(new MixedGradGradIntegrator);
    a01.AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(coeff_r_inv_hat));
    a01.Assemble();
    a01.EliminateTrialDofs(ess_bdr_1, stream, b0);
    a01.EliminateTestDofs(ess_bdr_0);    
    a01.Finalize();
    A01 = a01.ParallelAssemble();

    b0.ParallelAssemble(B0);
}

//Boundary condition for stream
double boundary_stream_f(const Vector &x){
    double x_2 = pow(x(0), 2);
    if (x(0) < 1.)
        return x_2*(2-x_2);
    return 1.;
}

//Boundary condition for stream at inlet
double boundary_stream_in_f(const Vector &x){
    double x_2 = pow(x(0), 2);
    return x_2*(2-x_2);
}
