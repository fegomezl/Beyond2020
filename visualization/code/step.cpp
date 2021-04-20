#include "header.h"

void Artic_sea::time_step(){
    //Update time
    t += config.dt;

    //Check for last iteration
    last = (t >= config.t_final - 1e-8*config.dt);

    //Update solution
    function.SetTime(t);
    x->ProjectCoefficient(function);

    //Output from the solution
    paraview_out->SetCycle(iteration);
    paraview_out->SetTime(t);
    paraview_out->Save();
    
    //Print the system state 
    double percentage = 100*(t - config.t_init)/(config.t_final - config.t_init);
    string progress = to_string((int)percentage)+"%";
    if (config.master){    
        cout.precision(4);
        cout << left << setw(12)
             << iteration << setw(12)
             << t << setw(12)
             << progress << "\r";
        cout.flush();
    }
}
