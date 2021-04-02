#include "header.h"

void Artic_sea::output_results(){
    //Print general information of the program
    if (master){
        cout << "\nSize: " << size << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << max(refinements-serial_refinements,0) << "\n"
             << "Total refinements: " << refinements << "\n";
    }

    //Output to Paraview
    ParaViewDataCollection paraview_out("graph", pmesh);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.SetCycle(0);
    paraview_out.SetTime(0.0);
    paraview_out.SetLevelsOfDetail(order);
    paraview_out.RegisterField("Temperature", x);
    paraview_out.Save();
}
