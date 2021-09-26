#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL),
    fec(NULL), fec_v(NULL),
    fespace(NULL), fespace_v(NULL),
    w(NULL), psi(NULL), v(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    solve_system();
    total_time = toc();
    output_results();
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    delete fec;
    delete fec_v;
    delete fespace;
    delete fespace_v;
    delete w;
    delete psi;
    delete v;
    delete M;
    delete D;
    delete C;
    delete Ct;
    if (config.master) cout << "Memory deleted!\n";
}
