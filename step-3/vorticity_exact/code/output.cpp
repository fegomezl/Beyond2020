#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nScalar size: " << size << "\n"
             << "Vectorial size: " << size_v << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "W total mean error: " << actual_error_w/(height*(out_rad-int_rad)) << "\n"
             << "Psi total mean error: " << actual_error_psi/(height*(out_rad-int_rad)) << "\n";
    }

    //Print visual results to Paraview
    if (config.last){
        ParaViewDataCollection paraview_out("results/graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.SetLevelsOfDetail(config.order);
        paraview_out.SetHighOrderOutput(true);
        paraview_out.RegisterField("StreamFunction(r)", psi);
        paraview_out.RegisterField("Vorticity(r)", w);
        paraview_out.RegisterField("Velocity(r)", v);
        paraview_out.Save();
    }
}
