#include "header.h"

//Initial condition

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_H1, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_H1, BlockVector &X):
    config(config),
    TimeDependentOperator(2*fespace_H1.GetTrueVSize(), 0.),
    fespace_H1(fespace_H1), fespace_ND(fespace_ND),
    block_offsets_H1(block_offsets_H1),
    m0(NULL), m1(NULL), 
    k0(NULL), k1(NULL), 
    b0(NULL), b1(NULL), 
    M0(NULL), M1(NULL), 
    M0_o(NULL), M1_o(NULL), 
    M0_e(NULL), M1_e(NULL), 
    K0(NULL), K1(NULL), 
    T0(NULL), T1(NULL),
    T0_e(NULL), T1_e(NULL),
    B0(NULL), B1(NULL), 
    B0_dt(NULL), B1_dt(NULL),
    Z0(&fespace_H1), Z1(&fespace_H1),
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
    Array<int> ess_bdr_0(attributes);
    ess_bdr_0 = 0;
    ess_bdr_0[0] = 0;  ess_bdr_0[1] = 0;
    ess_bdr_0[2] = 0;  ess_bdr_0[3] = 0;
    ess_bdr_0[4] = 1;  ess_bdr_0[5] = 0;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);

    Array<int> ess_bdr_1(attributes);
    ess_bdr_1 = 0;
    ess_bdr_1[0] = 0;  ess_bdr_1[1] = 0;
    ess_bdr_1[2] = 0;  ess_bdr_1[3] = 0;
    ess_bdr_1[4] = 1;  ess_bdr_1[5] = 0;
    fespace_H1.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    //Define initial condition
    ParGridFunction enthalphy(&fespace_H1);
    ConstantCoefficient coeff_enthalphy_0(1.1);
    enthalphy.ProjectCoefficient(coeff_enthalphy_0);
    enthalphy.GetTrueDofs(X.GetBlock(0));

    ParGridFunction salinity(&fespace_H1);
    ConstantCoefficient coeff_salinity_0(0.);
    salinity.ProjectCoefficient(coeff_salinity_0);
    salinity.GetTrueDofs(X.GetBlock(1));

    //Create mass matrix
    m0 = new ParBilinearForm(&fespace_H1);
    m0->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m0->Assemble();
    m0->Finalize();
    M0 = m0->ParallelAssemble();
    M0_e = M0->EliminateRowsCols(ess_tdof_0);
    M0_o = m0->ParallelAssemble();

    m1 = new ParBilinearForm(&fespace_H1);
    m1->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m1->Assemble();
    m1->Finalize();
    M1 = m1->ParallelAssemble();
    M1_e = M1->EliminateRowsCols(ess_tdof_1);
    M1_o = m1->ParallelAssemble();

    //Configure M solver
    M0_prec.SetPrintLevel(0);
    M0_prec.SetOperator(*M0);
    M0_solver.SetTol(config.reltol_conduction);
    M0_solver.SetAbsTol(config.abstol_conduction);
    M0_solver.SetMaxIter(config.iter_conduction);
    M0_solver.SetPrintLevel(0);
    M0_solver.SetPreconditioner(M0_prec);
    M0_solver.SetOperator(*M0);

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
