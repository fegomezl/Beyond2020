#include "header.h"

//Rotational functions
double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

//Physical properties (in H,S)
double FusionPoint(const double S);
double Phase(const double H);
double HeatDiffusivity(const double H);
double SaltDiffusivity(const double H);
double FusionDiffusivity(const double H, const double S);
double InversePermeability(const double H); 
double ExpansivityEnthalphy(const double H, const double S);
double ExpansivitySalinity(const double H, const double S);
double Buoyancy(const double H, const double S);

//Relationship between variables
double TStoHS(const double T, const double S);
double HStoTS(const double H, const double S);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    enthalphy = new ParGridFunction(fespace_H1);
    salinity = new ParGridFunction(fespace_H1);
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
    flow_oper = new Flow_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1);

    //Solve the initial velocity field
    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Set initial global state
    enthalphy->Distribute(X.GetBlock(0));
    salinity->Distribute(X.GetBlock(1));
    vorticity->Distribute(Y.GetBlock(0));
    stream->Distribute(Y.GetBlock(1));
    velocity->Distribute(Velocity);
    r_velocity->Distribute(r_Velocity);

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
    paraview_out->RegisterField("Enthalphy", enthalphy);
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

double Phase(const double H){
    return 0.25*(H+1)*tanh(5*JumpScale*(H+1))
         - 0.25*(H-1)*tanh(5*JumpScale*(H-1))
         + 0.5;
}

double HeatDiffusivity(const double H){
    double Dimensions = pow(LenghtScale, 2)/TimeScale;
    double liquid = 8.36;   //k_l/c_l
    double solid  = 57.63;  //k_s/c_s
    return (solid + (liquid-solid)*Phase(H))*Dimensions;
    return solid *(0.5*(1 - tanh(5*JumpScale*(H+1))))
         + liquid*(0.5*(1 + tanh(5*JumpScale*(H-1))));
}

double SaltDiffusivity(const double H){
    double Dimensions = pow(LenghtScale, 2)/TimeScale;
    double liquid = 0.1;    //d_l
    double solid  = 0.;     //d_s
    return (solid + (liquid-solid)*Phase(H))*Dimensions;
}

double FusionDiffusivity(const double H, const double S){
    double Dimensions = pow(LenghtScale, 2)/TimeScale;
    double a = 0.6037;
    double b = 0.00058123;
    double liquid = 0.2000;    //k_l/(L/2)
    double solid  = 0.6718;    //k_s/(L/2)
    return (a + 3*b*pow(S, 2))*(solid + (liquid-solid)*Phase(H))*Dimensions;
}

double InversePermeability(const double H){
    return Epsilon + pow(1-Phase(H), 2)/(pow(Phase(H), 3) + Epsilon);
} 

double ExpansivityEnthalphy(const double H, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = -562.4,   b0 = 45.5,
           a1 = 885.4,    b1 = -64.14,
           a2 = -598.207,
           a3 = 221.1715;
    return Phase(H)*((a0 + a1*(H-1) + a2*pow(H-1, 2) + a3*pow(H-1, 3))*S 
                   + (b0 + b1*(H-1))*pow(abs(S), 1.5))*Scale;
} 

double ExpansivitySalinity(const double H, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = 2697.0,   b0 = -89.7,  c0 = 32,
           a1 = -562.4,   b1 = 68.3,
           a2 = 442.7,    b2 = -48.10,
           a3 = -199.402,
           a4 = 55.29286;
    return Phase(H)*((a0 + a1*(H-1) + a2*pow(H-1, 2) + a3*pow(H-1, 3) + a4*pow(H-1, 4)) +
           (b0 + b1*(H-1) + b2*pow(H-1, 2))*pow(abs(S), 0.5) + (c0)*S)*Scale;
}

double Buoyancy(const double H, const double S){
    double Scale = pow(LenghtScale*TimeScale, -1);
    double a0 = 2697.0,   b0 = -59.8,    c0 = 16,
           a1 = -562.4,   b1 = 45.5,
           a2 = 442.7,    b2 = -32.1,
           a3 = -199.402,
           a4 = 55.29286;
    return Phase(H)*((a0 + a1*(H-1) + a2*pow(H-1, 2) + a3*pow(H-1, 3) + a4*pow(H-1, 4))*S +
           (b0 + b1*(H-1) + b2*pow(H-1, 2))*pow(abs(S), 1.5) + (c0)*pow(S, 2))*Scale;

}

double TStoHS(const double T, const double S){
    double T_rel = T - FusionPoint(S);
    double solid  = 0.0116;     //c_s/(L/2)
    double liquid = 0.0239;     //c_l/(L/2)
    return (solid *T_rel-1)*0.5*(1-tanh(5*JumpScale*T_rel))
         + (liquid*T_rel+1)*0.5*(1+tanh(5*JumpScale*T_rel));
}

double HStoTS(const double H, const double S){
    double solid  = 86.2069;    //(L/2)/c_s
    double liquid = 41.8410;    //(L/2)/c_l
    double T_rel = (solid *(H+1))*0.5*(1-tanh(5*JumpScale*(H+1)))
                 + (liquid*(H-1))*0.5*(1+tanh(5*JumpScale*(H-1)));
    return T_rel + FusionPoint(S);
} 
