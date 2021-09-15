#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    u(exact),
    r(rf),
    pmesh(NULL),
    fec(NULL),
    fespace(NULL),
    a(NULL),
    b(NULL),
    x(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    solve_system();
    output_results();
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    delete fec;
    delete fespace;
    delete a;
    delete b;
    delete x;
    if (config.master) cout << "Memory deleted!\n";
}
