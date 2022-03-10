#include "header.h"

//Initial conditions
double initial_temperature_f(const Vector &x);
double initial_salinity_f(const Vector &x);

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace_L2, ParFiniteElementSpace &fespace_ND, int dim, int attributes, Array<int> block_offsets_L2, BlockVector &X):
    config(config),
    TimeDependentOperator(2*fespace_L2.GetTrueVSize(), config.t_init),
    fespace_L2(fespace_L2),
    block_offsets_L2(block_offsets_L2),
    ess_bdr_0(attributes), ess_bdr_1(attributes),
    temperature(&fespace_L2), salinity(&fespace_L2),
    phase(&fespace_L2), 
    heat_inertia(&fespace_L2), heat_diffusivity(&fespace_L2), salt_diffusivity(&fespace_L2), 
    rvelocity(&fespace_ND), 
    coeff_r(r_f), 
    coeff_zero(dim, zero_f),
    coeff_rM(coeff_r, coeff_r),
    coeff_rD0(coeff_r, coeff_r),
    coeff_rD1(coeff_r, coeff_r),
    coeff_rV(&rvelocity), coeff_rMV(coeff_r, coeff_zero), 
    coeff_dPdT(coeff_zero, coeff_zero), coeff_dT_2(coeff_zero, coeff_zero),
    M0(NULL), M1(NULL), 
    K0(NULL), K1(NULL), 
    T0(NULL), T1(NULL),
    B0(NULL), B1(NULL), 
    B0_dt(&fespace_L2), B1_dt(&fespace_L2),
    Z0(&fespace_L2), Z1(&fespace_L2),
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
    //                            0

    ess_bdr_0 = 0;
    ess_bdr_0 [0] = 0;   ess_bdr_0 [1] = 0;   
    ess_bdr_0 [2] = 0;   ess_bdr_0 [3] = 0;
    ess_bdr_0 [4] = 1;   ess_bdr_0 [5] = 0;
    fespace_L2.GetEssentialTrueDofs(ess_bdr_0, ess_tdof_0);

    ess_bdr_1 = 0;
    ess_bdr_1 [0] = 0;     ess_bdr_1 [1] = 0;
    ess_bdr_1 [2] = 0;     ess_bdr_1 [3] = 0;
    ess_bdr_1 [4] = 1;     ess_bdr_1 [5] = 0;
    fespace_L2.GetEssentialTrueDofs(ess_bdr_1, ess_tdof_1);

    //Apply initial conditions
    if (!config.restart){
        FunctionCoefficient coeff_initial_temperature(initial_temperature_f);
        ConstantCoefficient coeff_boundary_temperature(InflowTemperature);
        temperature.ProjectCoefficient(coeff_initial_temperature);
        temperature.GetTrueDofs(X.GetBlock(0));
        
        FunctionCoefficient coeff_initial_salinity(initial_salinity_f);
        ConstantCoefficient coeff_boundary_salinity(InflowSalinity);
        salinity.ProjectCoefficient(coeff_initial_salinity);
        salinity.GetTrueDofs(X.GetBlock(1));
    }

    //Create mass matrix                                   
    ParBilinearForm m1(&fespace_L2);                 
    m1.AddDomainIntegrator(new MassIntegrator(coeff_r));  
    m1.Assemble();                                        
    m1.Finalize();                                        
    M1 = m1.ParallelAssemble();

    //Configure M solver
    M0_prec.SetPrintLevel(0);
    M0_solver.SetTol(config.reltol_conduction);
    M0_solver.SetAbsTol(config.abstol_conduction);
    M0_solver.SetMaxIter(config.iter_conduction);
    M0_solver.SetPrintLevel(0);

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
    bool NucleationRegion = (x(0) > RIn && x(1) > ZMax - NucleationHeight) || pow(x(0)-(RIn+ NucleationLength), 2) + pow(x(1)-(ZMax - NucleationHeight), 2) < pow(NucleationLength, 2);         //Wall with a small circle
    if (NucleationRegion)
        return NucleationTemperature;
    else
        return InitialTemperature;
}

double initial_salinity_f(const Vector &x){
  // if (x(0) < R_in && x(1) > Zmax - n_h/4)
  // return ((phi_out-phi_in)/n_h)*(x(1)-Zmax+n_h*1.2)+phi_in;

    bool NucleationRegion = (x(0) > RIn && x(1) > ZMax - NucleationHeight) || pow(x(0)-(RIn+ NucleationLength), 2) + pow(x(1)-(ZMax - NucleationHeight), 2) < pow(NucleationLength, 2);         //Wall with a small circle      
    if (NucleationRegion)
        return NucleationSalinity;
    else
        return InitialSalinity;
}
