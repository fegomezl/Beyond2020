#include "header.h"

Artic_sea::Artic_sea(Config config):config(config){}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    for (iteration; last; iteration++)
        time_step();
    output_results();
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    if (delete_fec) delete fec;
    delete fespace;
    delete x;
    delete oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted!\n";
}
