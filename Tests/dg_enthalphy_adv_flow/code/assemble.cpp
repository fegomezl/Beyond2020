#include "header.h"

//Rotational functions
double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

//Relationship between variables
double TtoH(const double T);
double HtoT(const double H);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    enthalphy = new ParGridFunction(fespace_L2);
    vorticity = new ParGridFunction(fespace_H1);
    stream = new ParGridFunction(fespace_H1);
    velocity = new ParGridFunction(fespace_ND);
    r_velocity = new ParGridFunction(fespace_ND);

    X = new HypreParVector(fespace_L2);
    Velocity = new HypreParVector(fespace_ND);
    r_Velocity = new HypreParVector(fespace_ND);

    //Integration setup
    int order_quad = max(2, 2*config.order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    //Initialize operators
    transport_oper = new Transport_Operator(config, *fespace_L2, *fespace_ND, dim, pmesh->bdr_attributes.Max(), *X);


    MPI_Barrier(MPI_COMM_WORLD);
    cout << "a\n";
    MPI_Barrier(MPI_COMM_WORLD);
    flow_oper = new Flow_Operator(config, *fespace_H1, *fespace_L2, *fespace_ND, dim, pmesh->bdr_attributes.Max(), block_offsets_H1);

    //Solve the initial velocity field
    flow_oper->SetParameters(*X);
    flow_oper->Solve(Y, *Velocity, *r_Velocity);

    //Set initial global state
    enthalphy->Distribute(X);
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

double TtoH(const double T){
    return 0.01906*T + (0.00731*T+1)*tanh(5*JumpScale*T);
}

double HtoT(const double H){
    return 61.51053*H+23.59584 
         + 18.95734*(H-1)*tanh(5*JumpScale*(H-1)) 
         + 42.55319*(H+1)*tanh(5*JumpScale*(H+1)); 
}

double Phase(const double H){
    return 0.5
         + 0.25*(H+1)*tanh(5*JumpScale*(H+1))
         + 0.25*(H-1)*tanh(5*JumpScale*(H-1));
}

double Density(const double H){
    if (H < 1)
        return 0.;
    else {
        double T = HtoT(H);
        return 2*T*(8-T);
    }
}
