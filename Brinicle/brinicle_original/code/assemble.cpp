#include "header.h"

//Usefull position functions
double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

//Physical properties (in T,S)
double FusionPoint(const double S);
double Phase(const double T, const double S);
double HeatInertia(const double T, const double S);
double HeatDiffusivity(const double T, const double S);
double SaltDiffusivity(const double T, const double S);
double Impermeability(const double T, const double S);
double Density(const double T, const double S);

//Initialize the solvers and the variables of the program
void Artic_sea::assemble_system(){

    if (!config.restart){
        //Define the temperature and salinity fields
        temperature = new ParGridFunction(fespace_H1);
        salinity = new ParGridFunction(fespace_H1);
    } else {
        //Read input temperature and salinity fields
        std::ifstream in;
        std::ostringstream oss;
        oss << std::setw(10) << std::setfill('0') << config.pid;
        std::string n_temperature = "results/restart/temperature_"+oss.str()+".gf";
        std::string n_salinity = "results/restart/salinity_"+oss.str()+".gf";

        in.open(n_temperature.c_str(),std::ios::in);
        temperature = new ParGridFunction(pmesh, in);
        in.close();
        temperature->GetTrueDofs(X.GetBlock(0));

        in.open(n_salinity.c_str(),std::ios::in);
        salinity = new ParGridFunction(pmesh, in);
        in.close();
        salinity->GetTrueDofs(X.GetBlock(1));
    }

    //Define other fields
    phase = new ParGridFunction(fespace_H1);
    vorticity = new ParGridFunction(fespace_H1);
    stream = new ParGridFunction(fespace_H1);
    velocity = new ParGridFunction(fespace_ND);
    rvelocity = new ParGridFunction(fespace_ND);

    rVelocity = new HypreParVector(fespace_ND);
    Velocity = new HypreParVector(fespace_ND);

    //Initialize operators
    transport_oper = new Transport_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1, X);
    flow_oper = new Flow_Operator(config, *fespace_H1, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1, X);

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
        (*phase)(ii) = Phase((*temperature)(ii), (*salinity)(ii));

    //Adimentionalize stream function
    if (config.rescale){
        double stream_local_max = stream->Max(), stream_max;
        double stream_local_min = stream->Min(), stream_min;
        MPI_Allreduce(&stream_local_max, &stream_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&stream_local_min, &stream_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        for (int ii = 0; ii < stream->Size(); ii++)
                (*stream)(ii) = ((*stream)(ii)-stream_min)/(stream_max-stream_min);
    }

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
    paraview_out->SetTime(t);
    paraview_out->Save();

    //Start program check
    if (config.master){
        cout << left << setw(12)
             << "--------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Dt" << setw(12)
             << "Time" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------------------------\n";

        ofstream out;
        out.open("results/progress.txt", std::ios::trunc);
        out << left << setw(12) 
            << "Step" << setw(12)
            << "Dt" << setw(12)
            << "Time" << setw(12)
            << "Progress" << "\n";
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
double FusionPoint(const double S){
    return -(constants.FusionPoint_a*S + constants.FusionPoint_b*pow(S, 3));
}

//Phase indicator (1 for liquid and 0 for solid)
double Phase(const double T, const double S){
    return 0.5*(1+tanh(5*EpsilonInv*(T-FusionPoint(S))));
}

//Coefficient for the mass term in the temperature equation
double HeatInertia(const double T, const double S){
    return constants.TemperatureMass_s + (constants.TemperatureMass_l-constants.TemperatureMass_l)*Phase(T, S);
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
    return Phase(T, S)*(
          (constants.Density_a0  + 
           constants.Density_a1*(T)  + 
           constants.Density_a2*pow(T, 2)      + 
           constants.Density_a3*pow(T, 3)      +
           constants.Density_a4*pow(T, 4))*(S) + 
          (constants.Density_b0  + 
           constants.Density_b1*(T)  +
           constants.Density_b2*pow(T, 2))*pow(abs(S), 1.5)  +
          (constants.Density_c0)*pow(S, 2)
          );
}
