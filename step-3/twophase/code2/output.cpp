#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nSize(RT): " << size_rt << "\n"
             << "Size(L2): " << size_l2 << "\n"
             << "Total Size: " << size_rt + size_l2 << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n";

    //Print numerical results of the program
        ofstream output;
        output.precision(4);
        if (config.refinements != 0)
            output.open("results/convergence.txt", std::ios::app);
        else{
            output.open("results/convergence.txt", std::ios::trunc);
            output << left << setw(16) 
                   << "h" << setw(16) 
                   << "V error" << setw(16)
                   << "P error" << "\n";
        }
        output << left << setw(16) 
               << h_min << setw(16) 
               << v_error << setw(16)
               << p_error << "\n";
        output.close();
    }

    //Print visual results to Paraview
    if (config.last){
        ParaViewDataCollection paraview_out("results/graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.SetLevelsOfDetail(config.order);
        paraview_out.SetHighOrderOutput(true);
        paraview_out.RegisterField("Velocity", v);
        paraview_out.RegisterField("Pressure", p);
        paraview_out.Save();
    }
}
