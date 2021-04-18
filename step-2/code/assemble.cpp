#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    const double reltol = 1e-4, abstol = 1e-4;
    t = config.t_init;     
    dt = config.dt_init;
    last = false;
    boundary.SetTime(config.t_init);

    //Set boundary conditions
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(boundary);
    x->ProjectBdrCoefficient(boundary, ess_bdr);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(*fespace, X, ess_bdr);

    //Set the ODE solver type
    switch (config.ode_solver_type){
        // MFEM explicit methods
        case 1: ode_solver = new ForwardEulerSolver; break;
        case 2: ode_solver = new RK2Solver(0.5); break; // midpoint method
        case 3: ode_solver = new RK3SSPSolver; break;
        case 4: ode_solver = new RK4Solver; break;
        // MFEM implicit L-stable methods
        case 5: ode_solver = new BackwardEulerSolver; break;
        case 6: ode_solver = new SDIRK23Solver(2); break;
        case 7: ode_solver = new SDIRK33Solver; break;
        // CVODE
        case 8:
           cvode = new CVODESolver(MPI_COMM_WORLD, CV_ADAMS);
           cvode->Init(*oper);
           cvode->SetSStolerances(reltol, abstol);
           cvode->SetMaxStep(dt);
           ode_solver = cvode; break;
        case 9:
           cvode = new CVODESolver(MPI_COMM_WORLD, CV_BDF);
           cvode->Init(*oper);
           cvode->SetSStolerances(reltol, abstol);
           cvode->SetMaxStep(dt);
           ode_solver = cvode; break;
        // ARKODE
        case 10:
        case 11:
           arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::EXPLICIT);
           arkode->Init(*oper);
           arkode->SetSStolerances(reltol, abstol);
           arkode->SetMaxStep(dt);
           if (config.ode_solver_type == 11) { arkode->SetERKTableNum(FEHLBERG_13_7_8); }
           ode_solver = arkode; break;
        case 12:
           arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::IMPLICIT);
           arkode->Init(*oper);
           arkode->SetSStolerances(reltol, abstol);
           arkode->SetMaxStep(dt);
           ode_solver = arkode; break;
        default:
           cout << "Unknown ODE solver type: " << config.ode_solver_type << '\n'
                << "Setting ODE to BackwardEulerSolver.\n";
                 ode_solver = new BackwardEulerSolver; break;
    }

    // Initialize MFEM integrators, SUNDIALS integrators are initialized above
    if (config.ode_solver_type < 8)
        ode_solver->Init(*oper);

    // Since we want to update the diffusion coefficient after every time step,
    // we need to use the "one-step" mode of the SUNDIALS solvers.
    if (cvode) 
        cvode->SetStepMode(CV_ONE_STEP); 
    if (arkode)
        arkode->SetStepMode(ARK_ONE_STEP); 

    //Open the paraview output and print initial state
    paraview_out = new ParaViewDataCollection("graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);
    paraview_out->SetCycle(0);
    paraview_out->SetTime(0);
    paraview_out->Save();

    //Start program check
    if (config.master)
        cout << left << setw(12)
             << "--------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Time" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------\n";
}

double theta(double x, double alpha){
    return exp(-x/alpha)/sqrt(x) - sqrt(M_PI/alpha)*erfc(sqrt(x/alpha));
}

double exact(const Vector &x, double t){
    double eta = pow(x.Norml2(),2)/(4*(alpha_s + alpha_l)*t);
    if (eta > lambda)
        return T_i - (T_i - T_f)*theta(eta, alpha_l)/theta(lambda, alpha_l);
    else
        return T_f - (T_i - T_f)*(theta(eta, alpha_s) - theta(lambda, alpha_s));
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    fespace(fespace),
    m(NULL),
    k(NULL),
    T(NULL),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm()),
    z(height)
{
    const double rel_tol = 1e-8;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Construct M
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble(0);
    m->FormSystemMatrix(ess_tdof_list, M);

    //Configure M solver
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(M);
    M_prec.SetType(HypreSmoother::Jacobi);

    //Configure T solver
    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(X, ess_bdr);
}
