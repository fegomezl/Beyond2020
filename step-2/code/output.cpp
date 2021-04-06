#include "header.h"

void Artic_sea::output_results(){
    if (config.master){
    //Print general information of the program
        cout << "\nSize: " << size << "\n"
             << "Mesh side: " << h_min << "\n"
             << "Serial refinements: " << config.serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements-config.serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "Total iterations: " << iteration << "\n";
    }
}
