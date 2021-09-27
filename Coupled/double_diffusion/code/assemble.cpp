#include "header.h"

//Rotational functions
double r_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);

//Fusion temperature dependent of salinity
double T_fun(const double &salinity);

//Variation on solid heat capacity
double delta_c_s_fun(const double &temperature, const double &salinity);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    theta = new ParGridFunction(fespace);
    phi = new ParGridFunction(fespace);
    phase = new ParGridFunction(fespace);

    cond_oper = new Conduction_Operator(config, *fespace, dim, pmesh->bdr_attributes.Max(), block_true_offsets, X);

    //Recover initial conditions
    theta->Distribute(&(X.GetBlock(0)));
    phi->Distribute(&(X.GetBlock(1)));

    //Calculate phases
    for (int ii = 0; ii < phase->Size(); ii++){
        double T_f = config.T_f + T_fun((*phi)(ii));
        (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
    }

    //Set the ODE solver type
    arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::IMPLICIT);
    arkode->Init(*cond_oper);
    arkode->SetSStolerances(config.reltol_sundials, config.abstol_sundials);
    arkode->SetMaxStep(dt);
    arkode->SetStepMode(ARK_ONE_STEP);
    ode_solver = arkode;

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    paraview_out = new ParaViewDataCollection(folder, pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", theta);
    paraview_out->RegisterField("Salinity", phi);
    paraview_out->RegisterField("Phase", phase);
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

double r_f(const Vector &x){
    return x(0);
}

void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

double T_fun(const double &salinity){
    return -(0.6037*salinity + 0.00058123*pow(salinity, 3));
}

double delta_c_s_fun(const double &temperature, const double &salinity){
    return -0.00313*salinity -
           0.00704*temperature -
           0.0000783*salinity*temperature +
           16.9*salinity*pow(temperature, -2);
}
