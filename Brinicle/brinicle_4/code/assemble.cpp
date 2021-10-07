#include "header.h"

//Rotational functions
double r_f(const Vector &x);

//Fusion temperature dependent of salinity
double T_fun(const double &salinity){
    double a = 0.6037;
    double b = 0.00058123;
    return a*salinity + b*pow(salinity, 3);
}

double t_a(const Vector &x){
    double a00 = 10.27542, a01 = -0.83195,
           a10 = -0.38667, a11 = 0.02801,
           a20 = 0.00624,  
           a30 = -0.00006;
    double T = -x(0);
    double S = x(1);

    if (x(0) < T_fun(x(1))  && x(0) > 2 && x(1) < 22.5)
        return (a00 + a10*T + a20*pow(T, 2) + a30*pow(T, 3))*S +
               (a01 + a11*T)*pow(abs(S), 1.5);                                  
    else
        return 0.;
}

double t_b(const Vector &x){
    double a1 = -8.882737925,
           a2 = -9.438982566,
           a3 = 6.417890826;
    double T = -x(0);
    double S = x(1);

    if (x(0) < T_fun(x(1))  && x(0) > 2 && x(1) < 22.5)
        return a1 + a2*T + a3*S;               // k_t = g*b_t/mu
    else
        return 0.;
}

double s_a(const Vector &x){
    double b00 = -2061.94264, b01 = 68.54083, b02 = -24.43573,
           b10 = 10.27542,    b11 = -1.24793,
           b20 = -0.19333,    b21 = 0.02101,
           b30 = 0.00208,
           b40 = -0.00001;
    double T = -x(0);
    double S = x(1);

    if (x(0) < T_fun(x(1))  && x(0) > 2 && x(1) < 22.5)
        return (b00 + b10*T + b20*pow(T, 2) + b30*pow(T, 3) + b40*pow(T, 4)) +
               (b01 + b11*T + b21*pow(T, 2))*pow(abs(S), 0.5) +
               (b02)*S;                                                         
    else
        return 0.;
}

double s_b(const Vector &x){
    double b1 = -1975.119148,
           b2 = 7.468676765,
           b3 = -13.5844212;
    double T = -x(0);
    double S = x(1);

    if (x(0) < T_fun(x(1))  && x(0) > 2 && x(1) < 22.5)
        return b1 + b2*T + b3*S;             
    else
        return 0.;
}

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;
    actual_error = 0;
    total_error = 0;

    //Convergence setup
    int order_quad = max(2, 2*config.order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));

    ParGridFunction ta(fespace);
    FunctionCoefficient T_a(t_a);
    ta.ProjectCoefficient(T_a);

    ParGridFunction tb(fespace);
    FunctionCoefficient T_b(t_b);
    tb.ProjectCoefficient(T_b);

    ParGridFunction sa(fespace);
    FunctionCoefficient S_a(s_a);
    sa.ProjectCoefficient(S_a);

    ParGridFunction sb(fespace);
    FunctionCoefficient S_b(s_b);
    sb.ProjectCoefficient(S_b);

    double t_error = tb.ComputeL2Error(T_a, irs)/ComputeGlobalLpNorm(2, T_a, *pmesh, irs); 
    double s_error = sb.ComputeL2Error(S_a, irs)/ComputeGlobalLpNorm(2, S_a, *pmesh, irs); 

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    paraview_out = new ParaViewDataCollection(folder, pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("T_a", &ta);
    paraview_out->RegisterField("T_b", &tb);
    paraview_out->RegisterField("S_a", &sa);
    paraview_out->RegisterField("S_b", &sb);
    paraview_out->SetCycle(vis_impressions);
    paraview_out->SetTime(t);
    paraview_out->Save();

    //Start program check
    if (config.master)
        cout << left << setw(12)
             << "----------------------------------------------------------------\n"
             << left << setw(12)
             << t_error << setw(12)
             << s_error << setw(12)
             << left << setw(12)
             << "\n----------------------------------------------------------------\n";

}

Conduction_Operator::Conduction_Operator(Config config, ParFiniteElementSpace &fespace, Array<int> ess_bdr):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    m(NULL), k(NULL), 
    M(NULL), M_e(NULL), M_0(NULL),
    K_0(NULL),
    T(NULL), T_e(NULL), 
    Z(&fespace),
    coeff_r(r_f), r_alpha(alpha, coeff_r),
    M_solver(fespace.GetComm()), T_solver(fespace.GetComm())
{

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();
    M_e = M->EliminateRowsCols(ess_tdof_list); 
    M_0 = m->ParallelAssemble();

    k = new ParBilinearForm(&fespace);
    k->AddDomainIntegrator(new DiffusionIntegrator(r_alpha));
    k->Assemble();
    k->Finalize();
    K_0 = k->ParallelAssemble();    

    //Configure M solver
    M_solver.SetTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_prec.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);

    M_prec.SetOperator(*M);
    M_solver.SetOperator(*M);

    //Configure T solver
    T_solver.SetTol(config.reltol_conduction);
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction);
    T_solver.SetPrintLevel(0);
    T_prec.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);
}

double r_f(const Vector &x){
    return x(0);
}
