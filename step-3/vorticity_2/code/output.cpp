#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nSize(Psi): " << size_psi << "\n"
             << "Size(W): " << size_w << "\n"
             << "Total Size: " << size_w + size_psi << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n";
    }

    //Print visual results to Paraview
    if (config.last){
        ParaViewDataCollection paraview_out("results/graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.SetLevelsOfDetail(config.order);
        paraview_out.SetHighOrderOutput(true);
        paraview_out.RegisterField("StreamFunction", psi);
        paraview_out.RegisterField("Vorticity", w);
        paraview_out.RegisterField("StreamGradient", v);
        paraview_out.Save();
    }
}
