#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;     
    dt = config.dt_init;
    last = false;
    actual_error = 0;
    total_error = 0;
    iterations_error = 0;

    //Set boundary conditions
    dir_bdr.SetSize(pmesh->bdr_attributes.Max());
    new_bdr.SetSize(pmesh->bdr_attributes.Max());
    dir_bdr = 0;     new_bdr = 0; 
    dir_bdr[2] = 1;  new_bdr[3] = 1;

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    ConstantCoefficient outer(0.);
    x->ProjectBdrCoefficient(outer, dir_bdr);
    x->ProjectCoefficient(outer);
    x->GetTrueDofs(X);

    //Define added boundary solution
    ParGridFunction boundary(fespace);
    ConstantCoefficient boundary_coeff(T_f);
    boundary.ProjectCoefficient(boundary_coeff);
    boundary.GetTrueDofs(Boundary);

    //Create operator
    oper = new Conduction_Operator(*fespace, X, dir_bdr, new_bdr, t);

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

double exact(const Vector &x, double t){
    double r_2 = pow(x(0),2) + pow(x(1),2);
    double R_2 = pow(int_rad, 2) + 355.2*T_i*int_rad*t;
    if (R_2 <= r_2)
        return T_f;
    else
        return T_f + 0.89*T_i*int_rad*log(R_2/r_2);
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, Array<int> dir_bdr, Array<int> new_bdr, double t_init):
    TimeDependentOperator(fespace.GetTrueVSize(), t_init),
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

    fespace.GetEssentialTrueDofs(dir_bdr, ess_tdof_list);

    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble(0);
    m->FormSystemMatrix(ess_tdof_list, M);

    ConstantCoefficient inner(T_i);
    f = new ParLinearForm(&fespace);
    f->AddBoundaryIntegrator(new DomainLFIntegrator(inner), new_bdr);
    f->Assemble();
    F = *(f->ParallelAssemble());

    //Configure M solver
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.iterative_mode = false;
    M_solver.SetRelTol(rel_tol);
    M_solver.SetAbsTol(0.);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
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
