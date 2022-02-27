#include "header.h"

//Initial condition
double temperature_0_f(const Vector &x);
double salinity_0_f(const Vector &x);

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_L2, BlockVector &X):
    config(config),
    TimeDependentOperator(2*fespace_L2.GetTrueVSize(), 0.),
    fespace_L2(fespace_L2), fespace_ND(fespace_ND),
    block_offsets_L2(block_offsets_L2),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    m0(NULL), m1(NULL), 
    k0(NULL), k1(NULL), 
    b0(NULL), b1(NULL), 
    M0(NULL), M1(NULL), 
    K0(NULL), K1(NULL), 
    T0(NULL), T1(NULL),
    B0(NULL), B1(NULL), 
    B0_dt(NULL), B1_dt(NULL),
    Z0(&fespace_L2), Z1(&fespace_L2),
    coeff_r(r_f), 
    coeff_zero(dim, zero_f),
    M0_solver(MPI_COMM_WORLD), M1_solver(MPI_COMM_WORLD), 
    T0_solver(MPI_COMM_WORLD), T1_solver(MPI_COMM_WORLD)
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
    //
    ess_bdr_0 = 0;
    ess_bdr_0[0] = 0;  ess_bdr_0[1] = 0;
    ess_bdr_0[2] = 0;  ess_bdr_0[3] = 0;
    ess_bdr_0[4] = 1;  ess_bdr_0[5] = 0;
    fespace_L2.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);

    ess_bdr_1 = 0;
    ess_bdr_1[0] = 0;  ess_bdr_1[1] = 0;
    ess_bdr_1[2] = 0;  ess_bdr_1[3] = 0;
    ess_bdr_1[4] = 1;  ess_bdr_1[5] = 0;
    fespace_L2.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    //Define initial condition
    ParGridFunction temperature(&fespace_L2);
    FunctionCoefficient coeff_temperature_0(temperature_0_f);
    ConstantCoefficient coeff_temperature_in(InflowTemperature);
    temperature.ProjectCoefficient(coeff_temperature_0);
    temperature.ProjectBdrCoefficient(coeff_temperature_0, ess_bdr_0);
    temperature.GetTrueDofs(X.GetBlock(0));

    ParGridFunction salinity(&fespace_L2);
    FunctionCoefficient coeff_salinity_0(salinity_0_f);
    ConstantCoefficient coeff_salinity_in(InflowSalinity);
    salinity.ProjectCoefficient(coeff_salinity_0);
    salinity.ProjectBdrCoefficient(coeff_salinity_0, ess_bdr_1);
    salinity.GetTrueDofs(X.GetBlock(1));

    //Create mass matrix (Mass equation)
    m1 = new ParBilinearForm(&fespace_L2);
    m1->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m1->Assemble();
    m1->Finalize();
    M1 = m1->ParallelAssemble();

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

double temperature_0_f(const Vector &x){
    if ((x(0) > R_in && x(1) > Z - NucleationHeight) || pow(x(0)-(R_in + NucleationLength), 2) + pow(x(1)-(Z - NucleationHeight), 2) < pow(NucleationLength, 2))
        return NucleationTemperature;
    else
        return InitialTemperature;
}

double salinity_0_f(const Vector &x){
    if ((x(0) > R_in && x(1) > Z - NucleationHeight) || pow(x(0)-(R_in + NucleationLength), 2) + pow(x(1)-(Z - NucleationHeight), 2) < pow(NucleationLength, 2))
        return NucleationSalinity;
    else
        return InitialSalinity;
}
