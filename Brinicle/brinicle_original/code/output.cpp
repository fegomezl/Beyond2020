#include "header.h"

//Print the final results
void Artic_sea::output_results(){

    //Save final state
    std::ofstream out;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << config.pid;
    std::string n_mesh = "results/restart/pmesh_"+oss.str()+".msh";
    std::string n_temperature = "results/restart/temperature_"+oss.str()+".gf";
    std::string n_salinity = "results/restart/salinity_"+oss.str()+".gf";

    out.open(n_mesh.c_str(),std::ios::out);
    pmesh->ParPrint(out);
    out.close();

    out.open(n_temperature.c_str(),std::ios::out);
    temperature->Save(out);
    out.close();

    out.open(n_salinity.c_str(),std::ios::out);
    salinity->Save(out);
    out.close();

    //Update the initial time for future simulations
    if (config.master){
        out.open("settings/parameters.txt", std::ios::app);
        out << t << "       #Initial_time\n";
        out.close();
    }

    //Print general information of the program
    if (config.master){
        cout << "\n\nSize (H1): " << size_H1 << "\n"
             << "Size (ND): " << size_ND << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "Total iterations: " << iteration << "\n"
             << "Total printing: " << vis_print << "\n"
             << "Total execution time: " << total_time << " s" << "\n";

        std::ofstream out;
        out.open("results/state.txt", std::ios::trunc);
        out << "Size (H1): " << size_H1 << "\n"
            << "Size (ND): " << size_ND << "\n"
            << "Mesh Size: " << h_min << "\n"
            << "Serial refinements: " << serial_refinements << "\n"
            << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
            << "Total refinements: " << config.refinements << "\n"
            << "Total iterations: " << iteration << "\n"
            << "Total printing: " << vis_print << "\n"
            << "Total execution time: " << total_time << " s" << "\n";
        out.close();
    }
}
