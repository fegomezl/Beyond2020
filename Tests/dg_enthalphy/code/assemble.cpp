#include "header.h"

//Auxiliar functions
double r_f(const Vector &x);
void zero_f(const Vector &x, Vector &f);

//Initial condition
double H_0(const Vector &x);

double H(const double T);
double T(const double H);

void Artic_sea::assemble_system(){
    //Initialize the system
    t = 0;
    dt = config.dt_init;
    last = false;
    vis_steps = config.vis_steps_max;
    vis_impressions = 0;

    //Define solution x
    x = new ParGridFunction(fespace);
    X = new HypreParVector(fespace);

    //Initialize operators
    oper = new Transport_Operator(config, *fespace, dim, pmesh->bdr_attributes.Max(), *X);

    //Read initial global state
    x->Distribute(X);

    //Set the ODE solver type
    arkode = new ARKStepSolver(MPI_COMM_WORLD, ARKStepSolver::IMPLICIT);
    arkode->Init(*oper);
    arkode->SetSStolerances(config.reltol_sundials, config.abstol_sundials);
    arkode->SetMaxStep(dt);
    arkode->SetStepMode(ARK_ONE_STEP);
    ode_solver = arkode;

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    paraview_out = new ParaViewDataCollection(folder, pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Enthalphy", x);
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

Transport_Operator::Transport_Operator(Config config, ParFiniteElementSpace &fespace, int dim, int attributes, Vector &X):
    TimeDependentOperator(fespace.GetTrueVSize(), 0.),
    config(config),
    fespace(fespace),
    ess_bdr(attributes), ess_bdr_l(attributes), ess_bdr_s(attributes),
    m(NULL), k(NULL), b(NULL), 
    M(NULL), K(NULL), T(NULL),
    B(NULL), B_dt(NULL),
    Z(&fespace),
    aux(&fespace),
    zero(dim, zero_f),
    coeff_r(r_f), 
    coeff_rA(coeff_r, coeff_r), coeff_rA_l(config.a_l, coeff_r), coeff_rA_s(config.a_s, coeff_r),
    M_solver(MPI_COMM_WORLD), T_solver(MPI_COMM_WORLD)
{
    //Set boundary conditions
    //
    //                  1
    //            /------------\
    //            |            |
    //    (0)    2|            |3
    //            |            |
    //            \------------/
    //                  0
    
    /*ess_bdr[0] = 1;  ess_bdr[1] = 1;
    ess_bdr[2] = 0;  ess_bdr[3] = 0;

    ess_bdr_l = 0; ess_bdr_l[1] = 1;

    ess_bdr_s = 0; ess_bdr_s[0] = 1;*/

    ess_bdr   = 0;
    ess_bdr_l = 0;
    ess_bdr_s = 0;

    //Define initial condition
    FunctionCoefficient H_0_coeff(H_0);
    aux.ProjectCoefficient(H_0_coeff);
    aux.GetTrueDofs(X);

    //Create mass matrix
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator(coeff_r));
    m->Assemble();
    m->Finalize();
    M = m->ParallelAssemble();

    //Configure M solver
    M_prec.SetPrintLevel(0);
    M_prec.SetOperator(*M);

    M_solver.SetTol(config.reltol_conduction);
    M_solver.SetAbsTol(config.abstol_conduction);
    M_solver.SetMaxIter(config.iter_conduction);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(*M);

    //Configure T solver
    T_solver.SetRelTol(config.reltol_conduction);
    T_solver.SetAbsTol(config.abstol_conduction);
    T_solver.SetMaxIter(config.iter_conduction);
    T_solver.SetKDim(10);
    T_solver.SetPrintLevel(0);
    T_prec.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);

    SetParameters(X);
}

double r_f(const Vector &x){
    return x(0);
}

void zero_f(const Vector &x, Vector &f){
    f(0) = 0.;
    f(1) = 0.;
}

double H_0(const Vector &x){
    return H_s + (H_l-H_s)*x(1)/Z;
}

double H(const double T){
    //Map option
    return 0.01906*T + (0.00731*T+1)*tanh(5*JumpScale*T);

    //Conditional option
    /*
    if (T >= 0)
        return 0.02637*T+1;
    else
        return 0.01175*T-1;
    */
}

double T(const double H){
    //Map option
    return 61.51053*H+23.59584 
         + 18.95734*(H-1)*tanh(5*JumpScale*(H-1)) 
         + 42.55319*(H+1)*tanh(5*JumpScale*(H+1)); 

    //Conditional option
    /*
    if (H >= 1)
        return 37.91469*(H-1);
    else if (H <= -1)
        return 85.10638*(H+1);
    else
        return 0.;
    */
}
