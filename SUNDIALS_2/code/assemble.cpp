#include "header.h"

void Artic_sea::assemble_system(){
    //Set boundary conditions
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  
    ConstantCoefficient boundary(T_f);

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    ConstantCoefficient x_0(T_i);
    x->ProjectCoefficient(x_0);
    x->ProjectBdrCoefficient(boundary, ess_bdr);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(fespace, X, pmesh->bdr_attributes.Max());//*****************************

    //Initialize the system
    const double reltol = 1e-4, abstol = 1e-4;//*******************************************************
    dt = config.dt_init;
    t = 0.;
    last = false;

    //*************************************************************************************************
    // 9. Define the ODE solver used for time integration.
    switch (config.ode_solver_type)
    {
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
    }
    // Initialize MFEM integrators, SUNDIALS integrators are initialized above
    if (config.ode_solver_type < 8) 
        ode_solver->Init(*oper);                 

    // Since we want to update the diffusion coefficient after every time step,
    // we need to use the "one-step" mode of the SUNDIALS solvers.
    if (cvode) { cvode->SetStepMode(CV_ONE_STEP); }
    if (arkode) { arkode->SetStepMode(ARK_ONE_STEP); }
    //************************************************************************** 

    //Open the paraview output and print initial state
    paraview_out = new ParaViewDataCollection("graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);

    //Start program check
    cout.precision(4);
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
    double r_2 = pow(x.Norml2(),2);
    if (r_2 > 4*lambda*(alpha_s + alpha_l)*t)
        return T_i - (T_i - T_f)*theta(r_2/(4*(alpha_s + alpha_l)*t),alpha_l)/theta(lambda,alpha_l);
    else
        return T_f - (T_i - T_f)*(theta(r_2/(4*(alpha_s + alpha_l)*t),alpha_s) - theta(lambda,alpha_s));
}

Conduction_Operator::Conduction_Operator(ParFiniteElementSpace *&fespace, const Vector &X, double b_size):
    fespace(fespace),
    TimeDependentOperator(fespace->GetTrueVSize(), 0.),
    M_solver(fespace->GetComm()),
    T_solver(fespace->GetComm()),
    z(height),
    m(NULL),
    k(NULL),
    T(NULL)
{
    const double rel_tol = 1e-8;

    Array<int> ess_bdr(b_size);
    ess_bdr = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    //Construct M
    m = new ParBilinearForm(fespace);
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

    SetParameters(X);
}