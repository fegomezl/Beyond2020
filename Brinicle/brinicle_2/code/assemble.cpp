#include "header.h"

//Rotational functions
double r_f(const Vector &x);
double inv_r(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);
void r_inv_hat_f(const Vector &x, Vector &f);

//Fusion temperature dependent of salinity
double T_fun(const double &salinity);

//Variation of parameters
double delta_c_s_fun(const double &temperature, const double &salinity);
double delta_k_s_fun(const double &temperature, const double &salinity);
double delta_l_s_fun(const double &temperature, const double &salinity);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = config.t_init;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    if (!config.restart){
        theta = new ParGridFunction(fespace);
        phi = new ParGridFunction(fespace);
    } else {
        std::ifstream in;
        std::ostringstream oss;
        oss << std::setw(10) << std::setfill('0') << config.pid;
        std::string n_theta = "results/restart/theta_"+oss.str()+".gf";
        std::string n_phi = "results/restart/phi_"+oss.str()+".gf";

        in.open(n_theta.c_str(),std::ios::in);
        theta = new ParGridFunction(pmesh, in);
        in.close();
        theta->GetTrueDofs(X.GetBlock(0));

        in.open(n_phi.c_str(),std::ios::in);
        phi = new ParGridFunction(pmesh, in);
        in.close();
        phi->GetTrueDofs(X.GetBlock(1));
    }

    w = new ParGridFunction(fespace);
    psi = new ParGridFunction(fespace);
    v = new ParGridFunction(fespace_v);
    rv = new ParGridFunction(fespace_v);
    phase = new ParGridFunction(fespace);

    rV = new HypreParVector(fespace_v);
    V = new HypreParVector(fespace_v);

    //Initialize operators
    cond_oper = new Conduction_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), block_true_offsets, X);
    flow_oper = new Flow_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), block_true_offsets, X);

    //Solve initial velocity field
    flow_oper->Solve(Z, *V, *rV);

    //Read initial global state
    theta->Distribute(&(X.GetBlock(0)));
    phi->Distribute(&(X.GetBlock(1)));
    w->Distribute(&(Z.GetBlock(0)));
    psi->Distribute(&(Z.GetBlock(1)));
    v->Distribute(V);
    rv->Distribute(rV);

    //Calculate phases
    for (int ii = 0; ii < phase->Size(); ii++){
        double T_f = config.T_f + T_fun((*phi)(ii));
        (*phase)(ii) = 0.5*(1 + tanh(5*config.invDeltaT*((*theta)(ii) - T_f)));
    }

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
    paraview_out->RegisterField("Salinity", phi);
    paraview_out->RegisterField("Phase", phase);
    paraview_out->RegisterField("Vorticity", w);
    paraview_out->RegisterField("Stream", psi);
    paraview_out->RegisterField("Velocity", v);
    paraview_out->RegisterField("rVelocity", rv);
    paraview_out->SetCycle(vis_impressions);
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

double r_f(const Vector &x){
    return x(0);
}

double inv_r(const Vector &x){
    return pow(x(0) + 1E-9, -1);
}

void rot_f(const Vector &x, DenseMatrix &f){
    f(0,0) = 0.;  f(0,1) = 1.;
    f(1,0) = -1.; f(1,1) = 0.;
}

void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

void r_inv_hat_f(const Vector &x, Vector &f){
    f(0) = pow(x(0), -1);
    f(1) = 0.;
}

double T_fun(const double &salinity){
    double a = 0.6037;
    double b = 0.00058123;
    return -(a*salinity + b*pow(salinity, 3));
}

double delta_c_s_fun(const double &temperature, const double &salinity){
    double a = 0.00692;
    double b = -0.00307;
    double c = 0.0000768;
    double d = 16.6;
    return a*temperature +
           b*salinity +
           c*temperature*salinity +
           d*salinity*pow(temperature, -2);
}

double delta_k_s_fun(const double &temperature, const double &salinity){
    double a = 36.7;
    return a*salinity/temperature;
}

double delta_l_s_fun(const double &temperature, const double &salinity){
    double a = -1.94;
    double b = 16.6;
    double c = -0.0035;
    return a*temperature +
           b*salinity +
           c*salinity/temperature;
}
