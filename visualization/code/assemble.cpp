#include "header.h"

void Artic_sea::assemble_system(){
    //Initialize the system
    t = config.t_init;     
    last = false;
    function.SetTime(t);

    //Define solution x and apply initial conditions
    x = new ParGridFunction(fespace);
    x->ProjectCoefficient(function);

    //Open the paraview output and print initial state
    paraview_out = new ParaViewDataCollection("graph", pmesh);
    paraview_out->SetDataFormat(VTKFormat::BINARY);
    paraview_out->SetLevelsOfDetail(config.order);
    paraview_out->RegisterField("Temperature", x);
    paraview_out->SetCycle(0);
    paraview_out->SetTime(t);
    paraview_out->Save();

    //Start program check
    if (config.master)
        cout << left << setw(12)
             << "--------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Time" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------------------------\n";
}

double theta(double x, double alpha){
    return exp(-x/alpha)/sqrt(x) - sqrt(M_PI/alpha)*erfc(sqrt(x/alpha));
}

double Function(const Vector &x, double t){
    double eta = pow(x.Norml2(),2)/(4*(alpha_s + alpha_l)*t);
    if (eta > lambda)
        return T_i - (T_i - T_f)*theta(eta, alpha_l)/theta(lambda, alpha_l);
    else
        return T_f - (T_i - T_f)*(theta(eta, alpha_s) - theta(lambda, alpha_s));
}
