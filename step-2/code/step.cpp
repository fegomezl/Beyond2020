#include "header.h"

void Artic_sea::time_step(){
    //Check for last iteration
    last = (t + dt >= config.t_final - dt/2.);

    //Perform the time_step
    ode_solver->Step(X, t, dt);

    //Print the system state, if neccesary
    if (last || (iteration % vis_steps) == 0){
        if (config.master) 
            cout << "step " << iteration 
                 << ", t = " << t << ".\n";
        x->SetFromTrueDofs(X);
        paraview_out->SetCycle(iteration);
        paraview_out->SetTime(t);
        paraview_out->Save();
    }

    //Update the operator
    oper->SetParameters(X);
}
