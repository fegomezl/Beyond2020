#include "header.h"

double initial_conditions(const Vector &x){
    return 0.1*(pow(x(0),2)+pow(x(1),2));
}

void Artic_sea::assemble_system(){

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    FunctionCoefficient x_0(initial_conditions);
    x->ProjectCoefficient(x_0);
    x->GetTrueDofs(X);

    //Create operator
    oper = new Conduction_Operator(*fespace, config.alpha, config.kappa, X);

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
