#include "header.h"

void Artic_sea::assemble_system(){

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    FunctionCoefficient x_0(initial_conditions);
    x->ProjectCoefficient(x_0);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(*fespace, config.t_init, config.alpha, config.kappa, X);

    //Set the ODE solver type
    switch (config.ode_solver_type){
        // Implicit L-stable methods
        case 1:  ode_solver = new BackwardEulerSolver; break;
        case 2:  ode_solver = new SDIRK23Solver(2); break;
        case 3:  ode_solver = new SDIRK33Solver; break;
        // Explicit methods
        case 11: ode_solver = new ForwardEulerSolver; break;
        case 12: ode_solver = new RK2Solver(0.5); break; // Midpoint method
        case 13: ode_solver = new RK3SSPSolver; break;
        case 14: ode_solver = new RK4Solver; break;
        case 15: ode_solver = new GeneralizedAlphaSolver(0.5); break;
        // Implicit A-stable methods (not L-stable)
        case 22: ode_solver = new ImplicitMidpointSolver; break;
        case 23: ode_solver = new SDIRK23Solver; break;
        case 24: ode_solver = new SDIRK34Solver; break;
        default:
           cout << "Unknown ODE solver type: " << config.ode_solver_type << '\n'
                << "Setting ODE to BackwardEulerSolver.\n";
                 ode_solver = new BackwardEulerSolver; break;
    }

    //Initialize the system
    ode_solver->Init(*oper);
    t = config.t_init;
    dt = config.dt_init;
    last = false;

    //Open the paraview output and print initial state
    x->SetFromTrueDofs(X);
    paraview_out = new ParaViewDataCollection("graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);
    paraview_out->SetCycle(0);
    paraview_out->SetTime(t);
    paraview_out->Save();
}

double initial_conditions(const Vector &X){
    return 0.1*(pow(X(0),2)+pow(X(1),2));
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, double t_init,
                                        double alpha, double kappa, const Vector &X):
    fespace(fespace),
    t_init(t_init),
    alpha(alpha),
    kappa(kappa),
    current_dt(0.),
    TimeDependentOperator(fespace.GetTrueVSize(), t_init),
    m(NULL),
    k(NULL),
    T(NULL),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm()),
    z(height)
{
    const double rel_tol = 1e-8;

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
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(M);

    //Configure T solver
    T_solver.iterative_mode = false;
    T_solver.SetRelTol(rel_tol);
    T_solver.SetAbsTol(0.);
    T_solver.SetMaxIter(100);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    //Update the bilinear operator K
    SetParameters(X);
}

void Conduction_Operator::SetParameters(const Vector &X){
    //Renovate k
    delete k;
    k = new ParBilinearForm(&fespace);

    //Arrange nonlinear operator
    ParGridFunction X_alpha(&fespace);
    X_alpha.SetFromTrueDofs(X);
    for (int ii = 0; ii < X_alpha.Size(); ii++)
        X_alpha(ii) = kappa + alpha*X_alpha(ii);
    GridFunctionCoefficient X_coeff(&X_alpha);

    //Create the new K
    k->AddDomainIntegrator(new DiffusionIntegrator(X_coeff));
    k->Assemble(0);
    k->FormSystemMatrix(ess_tdof_list, K);

    //Renovate T
    delete T;
    T = NULL;
}
