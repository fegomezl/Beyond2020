#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nSize (w): " << size_w << "\n"
             << "Size (psi): " << size_psi << "\n"
             << "Size (v): " << size_v << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n";

    //Print numerical results of the program
        /*ofstream output;
        output.precision(4);
        if (config.refinements != 0)
            output.open("results/convergence.txt", std::ios::app);
        else{
            output.open("results/convergence.txt", std::ios::trunc);
            output << left << setw(16) 
                   << "DOFs" << setw(16) 
                   << "h" << setw(16) 
                   << "L2 error" << "\n";
        }
        output << left << setw(16) 
               << size << setw(16) 
               << h_min << setw(16) 
               << l2_error << "\n";
        output.close();*/
    }

    //Print visual results to Paraview
    if (config.last){
        ParaViewDataCollection paraview_out("results/graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.SetLevelsOfDetail(config.order);
        paraview_out.RegisterField("Vorticity", x_w);
        paraview_out.RegisterField("StreamFunction", x_psi);
        paraview_out.RegisterField("StreamGradient", x_v);
        paraview_out.Save();
    }
}
