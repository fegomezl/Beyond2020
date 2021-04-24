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

double Function(const Vector &x, double t){
    double int_rad = 5;
    double r_2 = pow(x(0),2) + pow(x(1),2);
    double R_2 = pow(int_rad, 2) + 355.2*T_i*int_rad*t;
    if (R_2 <= r_2)
        return T_f;
    else
        return T_f + 0.89*T_i*int_rad*log(R_2/r_2);
}
