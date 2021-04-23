#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = config.t_init;     
    dt = config.dt_init;
    last = false;
    actual_error = 0;
    total_error = 0;
    iterations_error = 0;
    exact_coeff.SetTime(t);

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(exact_coeff);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(*fespace, X, t, pmesh->bdr_attributes.Max());

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

//Auxiliar functions
double theta(double x, double alpha){
    return exp(-x/alpha)/sqrt(x) - sqrt(M_PI/alpha)*erfc(sqrt(x/alpha));
}

double g(double x){
    if (x > lambda){
        double alpha_l_rel = alpha_l/(alpha_l + alpha_s);
        return alpha_l*exp(-x/alpha_l_rel)/(theta(x,alpha_l_rel)*sqrt(x));
    }
    else {
        double alpha_s_rel = alpha_s/(alpha_l + alpha_s);
        return alpha_s*exp(-x/alpha_s_rel)/sqrt(x);
    }
}

//Exact solution of the problem
double exact(const Vector &x, double t){
    double eta = pow(x.Norml2(),2)/(4*(alpha_s + alpha_l)*t);
    if (eta > lambda)
        return T_i - (T_i - T_f)*theta(eta, alpha_l)/theta(lambda, alpha_l);
    else
        return T_f - (T_i - T_f)*(theta(eta, alpha_s) - theta(lambda, alpha_s));
}

//Boundary conditions
double fnbc_1(const Vector &x, double t){
    double u = x(2);
    double r_2 = pow(x.Norml2(),2);
    double eta = r_2/(4*(alpha_s + alpha_l)*t);
    return (T_i - T_f)*(u/r_2)*g(eta);
}

double fnbc_2(const Vector &x, double t){
    double u = x(2);
    double r_2 = pow(x.Norml2(),2);
    double eta = r_2/(4*(alpha_s + alpha_l)*t);
    return -(T_i - T_f)*(u/r_2)*g(eta);
}

double fnbc_3(const Vector &x, double t){
    double u = sqrt(pow(x(0),2) + pow(x(1),2));
    double r_2 = pow(x.Norml2(),2);
    double eta = r_2/(4*(alpha_s + alpha_l)*t);
    return (T_i - T_f)*(u/r_2)*g(eta);
}

double fnbc_4(const Vector &x, double t){
    double u = sqrt(pow(x(0),2) + pow(x(1),2));
    double r_2 = pow(x.Norml2(),2);
    double eta = r_2/(4*(alpha_s + alpha_l)*t);
    return -(T_i - T_f)*(u/r_2)*g(eta);
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace &fespace, const Vector &X, double t_init, int bdr_size):
    TimeDependentOperator(fespace.GetTrueVSize(), t_init),
    fespace(fespace),
    nbc_1(fnbc_1),
    nbc_2(fnbc_2),
    nbc_3(fnbc_3),
    nbc_4(fnbc_4),
    m(NULL),
    k(NULL),
    f(NULL),
    T(NULL),
    M_solver(fespace.GetComm()),
    T_solver(fespace.GetComm()),
    z(height)
{
    const double rel_tol = 1e-8;

    //Set dirichlet boundary conditions
    ess_bdr.SetSize(bdr_size);
    ess_bdr = 0; 
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Set neumman boundary conditions
    nbc_marker_1.SetSize(bdr_size); nbc_marker_1 = 0;
    nbc_marker_2.SetSize(bdr_size); nbc_marker_2 = 0;
    nbc_marker_3.SetSize(bdr_size); nbc_marker_3 = 0;
    nbc_marker_4.SetSize(bdr_size); nbc_marker_4 = 0;
    nbc_marker_1[0] = 1;
    nbc_marker_2[1] = 1;
    nbc_marker_3[2] = 1;
    nbc_marker_4[3] = 1;

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

    SetParameters(X, t_init);
}
