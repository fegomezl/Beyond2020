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
    oper = new Conduction_Operator(fespace, X, pmesh->bdr_attributes.Max());

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
    dt = config.dt_init;
    t = 0;
    last = false;

    //Open the paraview output and print initial state
    paraview_out = new ParaViewDataCollection("graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);

    //Start program check
    if (config.master) 
        cout << "\n"
             << "---------------------------\n"
             << left << setw(8)
             << "Step" << setw(8)
             << "Time" << setw(8)
             << "Progress\n"
             << "---------------------------\n";
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
    current_dt(0.),
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
}
