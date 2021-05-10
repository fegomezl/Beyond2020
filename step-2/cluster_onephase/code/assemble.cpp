#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;
    actual_error = 0;
    total_error = 0;

    //Set boundary conditions
    ess_bdr.SetSize(pmesh->bdr_attributes.Max());
    ess_bdr = 1; ess_bdr[2]= 0; 

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(initial_f);
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

    if (config.master && config.print){
        ofstream output;
        string fname = "results/local_" + to_string(config.ode_solver_type) + "_" + to_string(config.order) + "_" + to_string(config.refinements) + ".txt";
        output.open(fname, ios::trunc);
        output.close();
    }
}

double rf(const Vector &x){
    return x(0);
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> ess_bdr):
  TimeDependentOperator(fespace.GetTrueVSize(), 0.),
  fespace(fespace),
  m(NULL),
  k(NULL),
  T(NULL),
  r(rf),
  M_solver(fespace.GetComm()),
  T_solver(fespace.GetComm()),
  z(&fespace)
{
  const double rel_tol = 1e-8;

  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  m = new ParBilinearForm(&fespace);
  m->AddDomainIntegrator(new MassIntegrator(r));
  m->Assemble(0);
  m->FormSystemMatrix(ess_tdof_list, M);

  //Configure M solver
  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
  M_solver.SetPreconditioner(M_prec);
  M_prec.SetType(HypreSmoother::Jacobi);
  M_solver.SetOperator(M);

  //Configure T solver
  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.);
  T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(0);
  T_solver.SetPreconditioner(T_prec);

  SetParameters(X);
}
