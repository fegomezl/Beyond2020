#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    r(rf),
    pmesh(NULL),
    fec_w(NULL), fec_psi(NULL), fec_v(NULL),
    fespace_w(NULL), fespace_psi(NULL), fespace_v(NULL),
    a_w(NULL), a_psi(NULL),
    b_w(NULL), b_psi(NULL),
    x_w(NULL), x_psi(NULL), x_v(NULL),
    B_w(NULL), B_psi(NULL),
    X_w(NULL), X_psi(NULL)
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
    delete fec_w, fec_psi, fec_v;
    delete fespace_w, fec_psi, fec_v;
    delete a_w, a_psi;
    delete b_w, b_psi;
    delete x_w, x_psi, x_v;
    delete B_w, B_psi;
    delete X_w, X_psi;
    if (config.master) cout << "Memory deleted!\n";
}
