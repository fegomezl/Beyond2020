#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    function(Function),
    pmesh(NULL), fec(NULL), fespace(NULL),
    x(NULL), paraview_out(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    for (iteration = 1; !last; iteration++)
        time_step();
    output_results();
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    delete fespace;
    delete fec;
    delete x;
    delete paraview_out;
    if (config.master) cout << "Memory deleted!\n";
}
