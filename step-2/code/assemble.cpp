#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = config.t_init;     
    dt = config.dt_init;
    last = false;
    actual_error = 0;
    total_error = 0;
    iterations_error = 0;
    boundary.SetTime(t);

    //Set boundary conditions
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(boundary);
    x->ProjectBdrCoefficient(boundary, ess_bdr);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(*fespace, X, ess_bdr, *x);

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
             << "--------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Time" << setw(12)
             << "Progress" << setw(12)
             << "Error"
             << left << setw(12)
             << "\n--------------------------------------------------\n";
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


Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr, ParGridFunction &x):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    fespace(fespace),
    m(NULL),
    k(NULL),
    f(NULL),
    T(NULL),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm()),
    z(height)
{
    const double rel_tol = 1e-8;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    ConstantCoefficient coeff(0.);
    f = new ParLinearForm(&fespace);
    f->AddDomainIntegrator(new DomainLFIntegrator(coeff));
    f->Assemble();

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

    SetParameters(X, ess_bdr, x);
}
