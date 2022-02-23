#include "header.h"

void Artic_sea::output_results(){
    //Print general information of the program
    if (config.master)
        cout << "\n\nSize(H1): " << size_H1 << "\n"
             << "Size(ND): " << size_ND << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "Total iterations: " << iteration << "\n"
             << "Total printing: " << vis_impressions << "\n"
             << "Total execution time: " << total_time << " s" << "\n";
}
