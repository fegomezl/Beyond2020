#include "header.h"

//Initial conditions
double initial_temperature_f(const Vector &x);
double initial_salinity_f(const Vector &x);

//Initialization of the solver
Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X):
    config(config),
    TimeDependentOperator(2*fespace_H1.GetTrueVSize(), 0.),
    fespace_H1(fespace_H1),
    block_offsets_H1(block_offsets_H1),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    temperature(&fespace_H1), salinity(&fespace_H1),
    relative_temperature(&fespace_H1), phase(&fespace_H1), 
    stream(&fespace_H1),
    heat_diffusivity(&fespace_H1), salt_diffusivity(&fespace_H1), 
    coeff_r(r_f), 
    coeff_zero(dim, zero_f),
    coeff_rot(dim, rot_f),
    M0(NULL), M1(NULL), 
    K0(NULL), K1(NULL), 
    T0(NULL), T1(NULL),
    Z0(&fespace_H1), Z1(&fespace_H1),
    M0_solver(MPI_COMM_WORLD), M1_solver(MPI_COMM_WORLD), 
    T0_solver(MPI_COMM_WORLD), T1_solver(MPI_COMM_WORLD)
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

    //Temperature boundary conditions
    ess_bdr_0 = 0;
    ess_bdr_0 [0] = 0;   ess_bdr_0 [1] = 0;   
    ess_bdr_0 [2] = 0;   ess_bdr_0 [3] = 0;
    ess_bdr_0 [4] = 1;   
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);

    //Salinity boundary conditions
    ess_bdr_1 = 0;
    ess_bdr_1 [0] = 0;     ess_bdr_1 [1] = 0;
    ess_bdr_1 [2] = 0;     ess_bdr_1 [3] = 0;
    ess_bdr_1 [4] = 1;     
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    //Apply initial conditions
    FunctionCoefficient coeff_initial_temperature(initial_temperature_f);
    ConstantCoefficient coeff_boundary_temperature(0.);
    temperature.ProjectCoefficient(coeff_initial_temperature);
    temperature.ProjectBdrCoefficient(coeff_boundary_temperature, ess_bdr_0);
    temperature.GetTrueDofs(X.GetBlock(0));
    
    FunctionCoefficient coeff_initial_salinity(initial_salinity_f);
    ConstantCoefficient coeff_boundary_salinity(0.);
    salinity.ProjectCoefficient(coeff_initial_salinity);
    salinity.ProjectBdrCoefficient(coeff_boundary_salinity, ess_bdr_1);
    salinity.GetTrueDofs(X.GetBlock(1));

    //Create mass matrix                                   
    ParBilinearForm m1(&fespace_H1);
    m1.AddDomainIntegrator(new MassIntegrator(coeff_r));
    m1.Assemble(); 
    m1.EliminateEssentialBC(ess_bdr_1, Operator::DIAG_ONE);
    m1.Finalize();                                        
    M1 = m1.ParallelAssemble();

    //Configure M solver
    M0_prec.SetPrintLevel(0);
    M0_solver.SetTol(config.reltol_conduction);
    M0_solver.SetAbsTol(config.abstol_conduction);
    M0_solver.SetMaxIter(config.iter_conduction);
    M0_solver.SetPrintLevel(0);
    M0_solver.SetPreconditioner(M0_prec);

    M1_prec.SetPrintLevel(0);
    M1_prec.SetOperator(*M1);
    M1_solver.SetTol(config.reltol_conduction);
    M1_solver.SetAbsTol(config.abstol_conduction);
    M1_solver.SetMaxIter(config.iter_conduction);
    M1_solver.SetPrintLevel(0);
    M1_solver.SetPreconditioner(M1_prec);
    M1_solver.SetOperator(*M1); 

    //Configure T solver 
    T0_prec.SetPrintLevel(0);
    T0_solver.SetTol(config.reltol_conduction);
    T0_solver.SetAbsTol(config.abstol_conduction);
    T0_solver.SetMaxIter(config.iter_conduction); 
    T0_solver.SetPrintLevel(0);                   
    T0_solver.SetPreconditioner(T0_prec);         
                                                  
    T1_prec.SetPrintLevel(0);                     
    T1_solver.SetTol(config.reltol_conduction);   
    T1_solver.SetAbsTol(config.abstol_conduction);
    T1_solver.SetMaxIter(config.iter_conduction); 
    T1_solver.SetPrintLevel(0);                   
    T1_solver.SetPreconditioner(T1_prec);
}

//Initial conditions
double initial_temperature_f(const Vector &x){
    bool NucleationRegion = (x(0) > 1. && x(1) > Z - NucleationHeight) || pow(x(0)-(1. + NucleationLength), 2) + pow(x(1)-(Z - NucleationHeight), 2) < pow(NucleationLength, 2);         //Wall with a small circle

    if (NucleationRegion)
        return NucleationTemperature;
    else if (x(1) > Z - NucleationHeight){
        return (Z-x(1))/NucleationHeight;
    } else
        return 1.;
}

double initial_salinity_f(const Vector &x){
    bool NucleationRegion = (x(0) > 1. && x(1) > Z - NucleationHeight) || pow(x(0)-(1. + NucleationLength), 2) + pow(x(1)-(Z - NucleationHeight), 2) < pow(NucleationLength, 2);         //Wall with a small circle

    if (NucleationRegion)
        return NucleationSalinity;
    else if (x(1) > Z - NucleationHeight){
        return (Z-x(1))/NucleationHeight; 
    } else
        return 1.;
}
