#include "header.h"

double r_f(const Vector &x);
double r_inv_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);
void r_inv_hat_f(const Vector &x, Vector &f);
void rot_f(const Vector &x, DenseMatrix &f);

double T_fun(const double &salinity);
double delta_c_s_fun(const double &temperature, const double &salinity);
double delta_k_s_fun(const double &temperature, const double &salinity);
double delta_l_s_fun(const double &temperature, const double &salinity);
double delta_rho_t_fun(const double &temperature, const double &salinity);
double delta_rho_p_fun(const double &temperature, const double &salinity);

void Artic_sea::assemble_system(){
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

    //Integration setup
    int order_quad = max(2, 2*config.order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    //Initialize operators
    cond_oper = new Conduction_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), block_true_offsets, X);
    flow_oper = new Flow_Operator(config, *fespace, *fespace_v, dim, pmesh->bdr_attributes.Max(), block_true_offsets, X);

    //Solve initial velocity field
    flow_oper->SetParameters(X);
    flow_oper->Solve(Z, *V, *rV);

    //Set initial state
    theta->Distribute(&(X.GetBlock(0)));
    phi->Distribute(&(X.GetBlock(1)));
    w->Distribute(&(Z.GetBlock(0)));
    psi->Distribute(&(Z.GetBlock(1)));
    v->Distribute(V);
    rv->Distribute(rV);
    
    //Calculate phases
    for (int ii = 0; ii < phase->Size(); ii++){
        double T_f = config.T_f + T_fun((*phi)(ii));
        (*phase)(ii) = 0.5*(1 + tanh(5*EpsilonInv*((*theta)(ii) - T_f)));
    }

    //Normalize stream
    if (config.rescale){
        double psi_local_max = psi->Max(), psi_max;
        double psi_local_min = psi->Min(), psi_min;
        MPI_Allreduce(&psi_local_max, &psi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&psi_local_min, &psi_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        for (int ii = 0; ii < psi->Size(); ii++)
                (*psi)(ii) = ((*psi)(ii)-psi_min)/(psi_max-psi_min);
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

double T_fun(const double &salinity){
    double a = 0.6037;
    double b = 0.00058123;
    return -(a*salinity + b*pow(salinity, 3));
}

double delta_rho_t_fun(const double &temperature, const double &salinity){
    //if (temperature < T_fun(salinity))
    //    return 0.;

    double a0 = -13.0,  b0 = 1.1,
           a1 = 0.5,    b1 = -0.04,
           a2 = -0.008,  
           a3 = 0.0001;
    return (a0 + a1*temperature + a2*pow(temperature, 2) + a3*pow(temperature, 3))*salinity +
           (b0 + b1*temperature)*pow(abs(salinity), 1.5);
}

double delta_rho_p_fun(const double &temperature, const double &salinity){
    //if (temperature < T_fun(salinity))
    //    return 0.;

    double a0 = 2617.9, b0 = -87.0, c0 = 31.0,
           a1 = -13.0,  b1 = 1.6,
           a2 = 0.2,    b2 = -0.03,
           a3 = -0.003,
           a4 = 0.00002;
    return (a0 + a1*temperature + a2*pow(temperature, 2) + a3*pow(temperature, 3) + a4*pow(temperature, 4)) +
           (b0 + b1*temperature + b2*pow(temperature, 2))*pow(abs(salinity), 0.5) +
           (c0)*salinity;
}



