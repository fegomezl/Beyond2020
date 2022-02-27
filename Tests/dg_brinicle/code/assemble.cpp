#include "header.h"

//Rotational functions
double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

//Physical properties (in H,S)
double FusionPoint(const double S);                           
double Phase(const double T, const double S);                 
double HeatInertia(const double T, const double S);       
double HeatDiffusivity(const double T, const double S);   
double SaltDiffusivity(const double T, const double S);       
double InversePermeability(const double T, const double S);   
double ExpansivityTemperature(const double T, const double S);
double ExpansivitySalinity(const double T, const double S);   

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    temperature = new ParGridFunction(fespace_L2);
    salinity = new ParGridFunction(fespace_L2);
    phase = new ParGridFunction(fespace_L2);
    vorticity = new ParGridFunction(fespace_H1);
    stream = new ParGridFunction(fespace_H1);
    velocity = new ParGridFunction(fespace_ND);
    r_velocity = new ParGridFunction(fespace_ND);

    Velocity = new HypreParVector(fespace_ND);
    r_Velocity = new HypreParVector(fespace_ND);

    //Integration setup
    int order_quad = max(2, 2*config.order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    //Initialize operators
    transport_oper = new Transport_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1, X);
    flow_oper = new Flow_Operator(config, *fespace_H1, *fespace_L2, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1);

    //Solve the initial velocity field
    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Set initial global state
    temperature->Distribute(X.GetBlock(0));
    salinity->Distribute(X.GetBlock(1));
    vorticity->Distribute(Y.GetBlock(0));
    stream->Distribute(Y.GetBlock(1));
    velocity->Distribute(Velocity);
    r_velocity->Distribute(r_Velocity);

    //Calculate phases
    for (int ii = 0; ii < phase->Size(); ii++)
        (*phase)(ii) = Phase((*temperature)(ii), (*salinity)(ii));

    //Set the ODE solver type
    arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::IMPLICIT);
    arkode->Init(*transport_oper);
    arkode->SetSStolerances(config.reltol_sundials, config.abstol_sundials);
    arkode->SetMaxStep(dt);
    arkode->SetStepMode(ARK_ONE_STEP);
    ode_solver = arkode;

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    paraview_out = new ParaViewDataCollection(folder, pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", temperature);
    paraview_out->RegisterField("Salinity", salinity);
    paraview_out->RegisterField("Vorticity", vorticity);
    paraview_out->RegisterField("Stream", stream);
    paraview_out->RegisterField("Velocity", velocity);
    paraview_out->RegisterField("r_Velocity", r_velocity);
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

double r_inv_f(const Vector &x){
    return pow(x(0), -1);
}

void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0), -1);
    f(1) = 0.;
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = 1.;
    f(1,0) = -1.; f(1,1) = 0.;
}

double FusionPoint(const double S){
    double a = 0.6037;
    double b = 0.00058123;
    return -(a*S + b*pow(S, 3));

}

double Phase(const double T, const double S){
    return 0.5*(1+tanh(5*JumpScale*(T-FusionPoint(S))));
}

double HeatInertia(const double T, const double S){
    double liquid = 0.0140;   //c_l/L
    double solid  = 0.0068;   //c_s/L
    return (solid + (liquid-solid)*Phase(T, S));
}

double HeatDiffusivity(const double T, const double S){ 
    double Scale = pow(LenghtScale, 2)/TimeScale;
    double liquid = 0.1 ;   //k_l/L
    double solid  = 0.3;   //k_s/L
    return (solid + (liquid-solid)*Phase(T, S))*Scale;
}

double SaltDiffusivity(const double T, const double S){
    double Scale = pow(LenghtScale, 2)/TimeScale;
    double liquid = 0.1;    //d_l
    double solid  = 0.;     //d_s
    return (solid + (liquid-solid)*Phase(T, S))*Scale;
}

double InversePermeability(const double T, const double S){
    return Epsilon + pow(1-Phase(T, S), 2)/(pow(Phase(T, S), 3) + Epsilon);
} 

double ExpansivityTemperature(const double T, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = -13.4,   b0 = 1.1,
           a1 = 0.5,     b1 = -0.04,
           a2 = -0.08,
           a3 = 0.0001;
    return Phase(T, S)*((a0 + a1*(T) + a2*pow(T, 2) + a3*pow(T, 3))*S 
                      + (b0 + b1*(T))*pow(abs(S), 1.5))*Scale;
} 

double ExpansivitySalinity(const double T, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = 2697.0,   b0 = -89.7,  c0 = 32.0,
           a1 = -13.4,    b1 = 1.6,
           a2 = 0.3,      b2 = -0.03,
           a3 = -0.003,
           a4 = 0.00002;
    return Phase(T, S)*((a0 + a1*(T) + a2*pow(T, 2) + a3*pow(T, 3) + a4*pow(T, 4)) +
                        (b0 + b1*(T) + b2*pow(T, 2))*pow(abs(S), 0.5) + (c0)*S)*Scale;
}

double Buoyancy(const double T, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = 2697.0,   b0 = -59.8,    c0 = 16.0,
           a1 = -13.4,    b1 = 1.1,
           a2 = 0.3,      b2 = -0.02,
           a3 = -0.003,
           a4 = 0.00002;
    return Phase(T, S)*((a0 + a1*(T) + a2*pow(T, 2) + a3*pow(T, 3) + a4*pow(T, 4))*S +
                        (b0 + b1*(T) + b2*pow(T, 2))*pow(abs(S), 1.5) + (c0)*pow(S, 2))*Scale;

}
