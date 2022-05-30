#include "header.h"

//Usefull position functions
double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

//Physical properties (in T,S)
double FusionPoint(const double T, const double S);
double Phase(const double T, const double S);
double HeatDiffusivity(const double T, const double S);
double SaltDiffusivity(const double T, const double S);
double Impermeability(const double T, const double S);
double Density(const double T, const double S);

//Initialize the solvers and the variables of the program
void Artic_sea::assemble_system(){

    //Define fields
    temperature = new ParGridFunction(fespace_H1);
    salinity = new ParGridFunction(fespace_H1);
    phase = new ParGridFunction(fespace_H1);
    vorticity = new ParGridFunction(fespace_H1);
    stream = new ParGridFunction(fespace_H1);
    velocity = new ParGridFunction(fespace_ND);
    rvelocity = new ParGridFunction(fespace_ND);

    rVelocity = new HypreParVector(fespace_ND);
    Velocity = new HypreParVector(fespace_ND);

    //Initialize operators
    transport_oper = new Transport_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1, X);
    flow_oper = new Flow_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1);

    //Solve initial velocity field
    flow_oper->SetParameters(X);
    flow_oper->Solve(Y, *Velocity, *rVelocity);

    //Set initial state
    temperature->Distribute(X.GetBlock(0));
    salinity->Distribute(X.GetBlock(1));
    vorticity->Distribute(Y.GetBlock(0));
    stream->Distribute(Y.GetBlock(1));
    velocity->Distribute(Velocity);
    rvelocity->Distribute(rVelocity);
    
    //Calculate phases
    for (int ii = 0; ii < phase->Size(); ii++)
        (*phase)(ii) = FusionPoint((*temperature)(ii), (*salinity)(ii));

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
    paraview_out->RegisterField("Phase", phase);
    paraview_out->RegisterField("Vorticity", vorticity);
    paraview_out->RegisterField("Stream", stream);
    paraview_out->RegisterField("Velocity", velocity);
    paraview_out->RegisterField("rVelocity", rvelocity);
    paraview_out->SetCycle(vis_print);
    paraview_out->SetTime(t*t_ref);
    paraview_out->Save();

    //Start program check
    if (config.master){
        cout.precision(4);
        cout << left << setw(12)
             << "--------------------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Printed" << setw(12)
             << "Dt(s)" << setw(12)
             << "Time(s)" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------------------------------------\n";
        cout << left << setw(12)
             << 0 << setw(12)
             << vis_print << setw(12)
             << dt*t_ref << setw(12)
             << t*t_ref  << setw(12)
             << "0%" << "\r";
    }
}

//Function for r
double r_f(const Vector &x){
    return x(0);
}

//Function for 1/r
double r_inv_f(const Vector &x){
    return pow(x(0), -1);
}

//Function for 0 (vector)
void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

//Function for (1/r)*r^ (r^ unitary vector)
void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0), -1);
    f(1) = 0.;
}

//Function for ( 0   1 )
//             (-1   0 )
void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = 1.;
    f(1,0) = -1.; f(1,1) = 0.;
}

//Fusion temperature at a given salinity
double FusionPoint(const double T, const double S){
    double T_ = T + ZeroTemperature;
    double S_ = S + ZeroSalinity;
    return (1-2*signbit(T_ref))*(T_ - (constants.FusionPoint_a*(S_) + constants.FusionPoint_b*pow(S_, 3)));
}

//Phase indicator (1 for liquid and 0 for solid)
double Phase(const double T, const double S){
    return 0.5*(1+tanh(5*EpsilonInv*FusionPoint(T, S)));
}

//Coefficient for the diffusion term in the temperature equation
double HeatDiffusivity(const double T, const double S){ 
    return constants.TemperatureDiffusion_s + (constants.TemperatureDiffusion_l-constants.TemperatureDiffusion_s)*Phase(T, S);
}

//Coefficient for the diffusion term in the salinity equation
double SaltDiffusivity(const double T, const double S){
    return constants.SalinityDiffusion_s + (constants.SalinityDiffusion_l-constants.SalinityDiffusion_s)*Phase(T, S);
}

//Inverse of the brinkman penalization permeability
double Impermeability(const double T, const double S){
    return Epsilon + pow(1-Phase(T, S), 2)/(pow(Phase(T, S), 3) + Epsilon);
} 

//Relative density of the fluid
double Density(const double T, const double S){
    double T_ = T + ZeroTemperature;
    double S_ = S + ZeroSalinity;
    return Phase(T, S)*(
          (constants.Density_a0                                +
           constants.Density_a1*(T_)                           +
           constants.Density_a2*pow(T_, 2)                     +
           constants.Density_a3*pow(T_, 3)                     +
           constants.Density_a4*pow(T_, 4))*(S_)               +
          (constants.Density_b0                                + 
           constants.Density_b1*(T_)                           +
           constants.Density_b2*pow(T_, 2))*pow(abs(S_), 1.5)  +
          (constants.Density_c0)*pow(S_, 2)
          );
}
