#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL),
    fec_psi(NULL), fec_w(NULL), fec_v(NULL),
    fespace_psi(NULL), fespace_w(NULL), fespace_v(NULL),
    block_offsets(3), block_true_offsets(3),
    f(NULL), g(NULL),
    m(NULL), d(NULL), c(NULL),
    M(NULL), D(NULL), C(NULL), Ct(NULL),
    A(NULL),
    w(NULL), psi(NULL), v(NULL)
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
    delete fec_psi, fec_w, fec_v;
    delete fespace_psi, fespace_w, fespace_v;
    delete f, g;
    delete m, c;
    delete M, C, Ct;  
    delete A;
    delete w, psi, v;
    if (config.master) cout << "Memory deleted!\n";
}
