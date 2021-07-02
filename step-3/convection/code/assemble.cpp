#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Set boundary conditions
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 0;  ess_bdr[0] = 1;  ess_bdr[1] = 1;

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(initial_f);
    x->ProjectBdrCoefficient(initial_f, ess_bdr);
    X = new HypreParVector(fespace);
    x->GetTrueDofs(*X);

    oper = new Conduction_Operator(*fespace, *X, ess_bdr);

    //Set the ODE solver type
    switch (config.ode_solver_type){
        // MFEM explicit methods
        case 1: ode_solver = new ForwardEulerSolver; break;
        case 2: ode_solver = new RK2Solver(0.5); break;
        case 3: ode_solver = new RK3SSPSolver; break;
        case 4: ode_solver = new RK4Solver; break;
        // MFEM implicit L-stable methods
        case 5: ode_solver = new BackwardEulerSolver; break;
        case 6: ode_solver = new SDIRK23Solver(2); break;
        case 7: ode_solver = new SDIRK33Solver; break;
        // CVODE
        case 8: cvode = new CVODESolver(MPI_COMM_WORLD, CV_ADAMS); break;
        case 9: cvode = new CVODESolver(MPI_COMM_WORLD, CV_BDF); break;
        // ARKODE
        case 10: arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::EXPLICIT); break;
        case 11: arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::EXPLICIT); break;
        case 12: arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::IMPLICIT); break;
        default:
                 cout << "Unknown ODE solver type: " << config.ode_solver_type << "\n"
                      << "Setting ODE to BackwardEulerSolver.\n";
                 config.ode_solver_type = 1;
                 ode_solver = new BackwardEulerSolver; break;
    }

    // Initialize ODE solver
    if (config.ode_solver_type < 8)
        ode_solver->Init(*oper);
    else if (cvode){
        cvode->Init(*oper);
        cvode->SetSStolerances(config.reltol, config.abstol);
        cvode->SetMaxStep(dt);
        cvode->SetStepMode(CV_ONE_STEP);
        ode_solver = cvode;
    }
    else if (arkode){
        arkode->Init(*oper);
        arkode->SetSStolerances(config.reltol, config.abstol);
        arkode->SetMaxStep(dt);
        arkode->SetStepMode(ARK_ONE_STEP);
        if (config.ode_solver_type == 11) arkode->SetERKTableNum(FEHLBERG_13_7_8);
        ode_solver = arkode;
    }

     //Open the paraview output and print initial state
    paraview_out = new ParaViewDataCollection("results/graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);
    paraview_out->SetCycle(vis_impressions);
    paraview_out->SetTime(t);
    paraview_out->Save();

    //Start program check
    if (config.master)
        cout << left << setw(12)
             << "--------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Dt" << setw(12)
             << "Time" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------------------------\n";

}

double initial(const Vector &x){
    double c = k_s*Zmax/(k_l + k_s);
    if (x(1) <= c)
        return -10*(1 - x(1)/c);
    else
        return 10*(x(1) - c)/(Zmax - c);

    //return (20*x(1)/Zmax - 10);
}

double rf(const Vector &x){
    return x(0);
}


double boundary_psi(const Vector &x){
    return x(0);
}

double f_rhs(const Vector &x){
    return (Rmax - Rmin)/2 - x(0);
}

double porous_constant(const Vector &x){
    double height = Zmax - Zmin;
    double mid_x = (Rmax + Rmin)/2;
    double mid_y = height/2;
    double sigma = (Rmax - Rmin)/5;

    double r_2 = pow(x(0) - mid_x, 2) + pow(x(1) - mid_y, 2);
    if (r_2 < pow(sigma, 2))
        return 1e+6;
    else
        return 0.1;
}


Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    fespace(fespace),
    m(NULL),
    k(NULL),
    t(NULL),
    aux(&fespace),
    aux_C(&fespace),
    aux_K(&fespace),
    r(rf),
    coeff_C(&aux_C), coeff_rC(r, coeff_C),
    coeff_K(&aux_K), coeff_rK(r, coeff_K), dt_coeff_rK(1., coeff_rK),
    coeff_L(&aux), coeff_rL(r, coeff_L),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm())
{
    const double rel_tol = 1e-8;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Configure M solver
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_prec.SetType(HypreSmoother::Jacobi);

    //Configure T solver
    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(X);
}

Flow_Operator::Flow_Operator(ParFiniteElementSpace &fespace_psi, ParFiniteElementSpace fespace_w,ParFiniteElementSpace fespace_v):
    fespace_psi(fespace_psi), fespace_w(fespace_w), fespace_v(fespace_v),
    block_offsets(3), block_true_offsets(3),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL), ct(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    A(NULL)
{
    //Create the block offsets
    block_offsets[0] = 0;
    block_offsets[1] = fespace_psi.GetVSize();
    block_offsets[2] = fespace_w.GetVSize();
    block_offsets.PartialSum();

    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace_psi.TrueVSize();
    block_true_offsets[2] = fespace_w.TrueVSize();
    block_true_offsets.PartialSum();

    //Initialize the corresponding vectors
    const char *device_config = "cpu";
    Device device(device_config);
    MemoryType mt = device.GetMemoryType();
    y.Update(block_offsets, mt); Y.Update(block_true_offsets, mt);
    b.Update(block_offsets, mt); B.Update(block_true_offsets, mt);
}