#include "header.h"

//Rotational functions
double r_f(const Vector &x);
double inv_r(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);
void r_inv_hat_f(const Vector &x, Vector &f);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define variables
    theta = new ParGridFunction(fespace);
    w = new ParGridFunction(fespace);
    psi = new ParGridFunction(fespace);
    v = new ParGridFunction(fespace_v);

    Theta = new HypreParVector(fespace);
    W = new HypreParVector(fespace);
    Psi = new HypreParVector(fespace);
    rV = new HypreParVector(fespace_v);
    V = new HypreParVector(fespace_v);

    //Initialize operators
    cond_oper = new Conduction_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), *Theta);
    flow_oper = new Flow_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), *Theta);

    //Solve initial velocity field
    flow_oper->Solve(*W, *Psi, *V, *rV);

    //Read initial global state
    theta->Distribute(Theta);
    w->Distribute(W);
    psi->Distribute(Psi);
    v->Distribute(V);

    //Normalize stream
    double psi_local_max = psi->Max(), psi_max;
    double psi_local_min = psi->Min(), psi_min;
    MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    for (int ii = 0; ii < psi->Size(); ii++)
            (*psi)(ii) = ((*psi)(ii)-psi_min)/(psi_max-psi_min);

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
    paraview_out->RegisterField("Vorticity", w);
    paraview_out->RegisterField("Stream", psi);
    paraview_out->RegisterField("Velocity", v);
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

double inv_r(const Vector &x){
    return pow(x(0), -1);
}

void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.; f(0,1) = -1.;
    f(1,0) = 1.; f(1,1) = 0.;
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0), -1);
    f(1) = 0.;
}
