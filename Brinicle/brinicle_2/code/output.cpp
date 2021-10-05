#include "header.h"

void Artic_sea::output_results(){
    //Save final state
    std::ofstream out;
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill('0') << config.pid;
    std::string n_mesh = "results/restart/pmesh_"+oss.str()+".msh";
    std::string n_theta = "results/restart/theta_"+oss.str()+".gf";
    std::string n_phi = "results/restart/phi_"+oss.str()+".gf";

    out.open(n_mesh.c_str(),std::ios::out);
    pmesh->ParPrint(out);
    out.close();

    out.open(n_theta.c_str(),std::ios::out);
    theta->Save(out);
    out.close();

    out.open(n_phi.c_str(),std::ios::out);
    phi->Save(out);
    out.close();

    if (config.master){
        out.open("settings/parameters.txt", std::ios::app);
        out << t << "       #Initial_time\n";
        out.close();
    }

    //Print general information of the program
    if (config.master)
        cout << "\n\nSize: " << 2*size << "\n"
             << "Vectorial size: " << size_v << "\n"
             << "Mesh Size: " << h_min << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n"
             << "Total iterations: " << iteration << "\n"
             << "Total printing: " << vis_impressions << "\n"
             << "Total execution time: " << total_time << " s" << "\n";
}
