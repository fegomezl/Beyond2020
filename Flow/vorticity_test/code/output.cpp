#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nScalar size: " << 2*size << "\n"
             << "Vectorial size: " << size_v << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "Total execution time: " << total_time << " s" << "\n";
    }

    //Print visual results to Paraview
    if (config.last){
        ParaViewDataCollection paraview_out("results/graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetLevelsOfDetail(config.order);
        paraview_out.RegisterField("Vorticity", w);
        paraview_out.RegisterField("Stream", psi);
        paraview_out.RegisterField("Velocity", v);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.Save();
    }
}
