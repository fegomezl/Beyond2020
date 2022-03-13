#include "header.h"

//Initial conditions
double initial_temperature_f(const Vector &x);
double initial_salinity_f(const Vector &x);

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, Vector &Temperature, Vector &Salinity):
    config(config),
    TimeDependentOperator(fespace_H1.GetTrueVSize(), config.t_init),
    fespace_H1(fespace_H1),
    block_offsets_H1(block_offsets_H1),
    ess_bdr(attributes),
    temperature(&fespace_H1), salinity(&fespace_H1),
    salt_diffusivity(&fespace_H1), 
    rvelocity(&fespace_ND), 
    coeff_r(r_f), 
    coeff_zero(dim, zero_f),
    coeff_rD(coeff_r, coeff_r),
    coeff_rV(&rvelocity),
    M(NULL),  
    M_e(NULL),
    M_o(NULL),
    K(NULL),  
    T(NULL),  
    T_e(NULL),
    B(NULL),  
    B_dt(&fespace_H1),
    Z(&fespace_H1),   
    M_solver(MPI_COMM_WORLD),
    T_solver(MPI_COMM_WORLD)
{
    //Define essential boundary conditions
    //   
    //                            1
    //            /-------------------------------\
    //            |                               |
    //            |                               |
    //            |                               | 
    //            |                               |
    //          2 |                               | 3
    //            |                               |
    //            |                               |
    //            |                               | 
    //            |                               |
    //            \-------------------------------/
    //                            0

    ess_bdr = 0;
    ess_bdr [0] = 0;   ess_bdr [1] = 0;   
    ess_bdr [2] = 0;   ess_bdr [3] = 0;
    fespace_H1.GetEssentialTrueDofs(ess_bdr, ess_tdof);

    //Apply initial conditions
    if (!config.restart){
        FunctionCoefficient coeff_initial_temperature(initial_temperature_f);
        temperature.ProjectCoefficient(coeff_initial_temperature);
        temperature.GetTrueDofs(Temperature);
        
        FunctionCoefficient coeff_initial_salinity(initial_salinity_f);
        ConstantCoefficient coeff_boundary_salinity(InflowSalinity);
        salinity.ProjectCoefficient(coeff_initial_salinity);
        salinity.ProjectBdrCoefficient(coeff_boundary_salinity, ess_bdr);
        salinity.GetTrueDofs(Salinity);
    }

    //Create mass matrix                                   
    ParBilinearForm m(&fespace_H1);                 
    m.AddDomainIntegrator(new MassIntegrator(coeff_r));  
    m.Assemble();                                        
    m.Finalize();                                        
    M = m.ParallelAssemble();
    M_e = M->EliminateRowsCols(ess_tdof);
    M_o = m.ParallelAssemble();

    //Set RHS
    ConstantCoefficient Zero(0.);
    ParLinearForm b(&fespace_H1);
    b.AddDomainIntegrator(new DomainLFIntegrator(Zero));
    b.Assemble();
    B = b.ParallelAssemble();

    //Configure M solver
    M_prec.SetPrintLevel(0);
    M_prec.SetOperator(*M);
    M_solver.SetTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(*M);                                                                                     

    //Configure T solver 
    T_prec.SetPrintLevel(0);                     
    T_solver.SetTol(config.reltol_conduction);   
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction); 
    T_solver.SetPrintLevel(0);                   
    T_solver.SetPreconditioner(T_prec);
}

//Initial conditions
double initial_temperature_f(const Vector &x){
    return -2.;
}

double initial_salinity_f(const Vector &x){
    if (x(0) > 0.8*RMax)
        return 1.;
    if (x(1) < ZMax/2)
        return 3.5;
    else
        return 6.5;
}
