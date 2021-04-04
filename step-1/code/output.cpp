#include "header.h"

void Artic_sea::output_results(){
    if (master){
    //Print general information of the program
        cout << "\nSize: " << size << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << max(refinements-serial_refinements,0) << "\n"
             << "Total refinements: " << refinements << "\n";

    //Print numerical results of the program
        ofstream output;
        output.precision(4);
        if (refinements != 0){
            output.open("data/convergence.txt", std::ios::app);
            output << left << setw(16) 
                   << size << setw(16) 
                   << h_min << setw(16) 
                   << l2_error << "\n";
        } else {
            output.open("data/convergence.txt", std::ios::trunc);
            output << left << setw(16) 
                   << "DOFs" << setw(16) 
                   << "h" << setw(16) 
                   << "L2 error" << "\n";
            output << left << setw(16) 
                   << size << setw(16) 
                   << h_min << setw(16) 
                   << l2_error << "\n";
        }
        output.close();
    }

    //Print visual results to Paraview
    if (last){
        ParaViewDataCollection paraview_out("graph", pmesh);
        paraview_out.SetDataFormat(VTKFormat::BINARY);
        paraview_out.SetCycle(0);
        paraview_out.SetTime(0.);
        paraview_out.SetLevelsOfDetail(order);
        paraview_out.RegisterField("Temperature", x);
        paraview_out.Save();
    }
}
