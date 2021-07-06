#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL),
    fec_w(NULL), fec_psi(NULL), fec_v(NULL),
    fespace_w(NULL), fespace_psi(NULL), fespace_v(NULL),
    block_offsets(3), block_true_offsets(3),
    g(NULL), f(NULL),
    m(NULL), d(NULL), c(NULL), ct(NULL),
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
    delete fec_w, fec_psi, fec_v;
    delete fespace_w, fespace_psi, fespace_v;
    delete g, f;
    delete m, c, ct;
    delete M, C, Ct;  
    delete A;
    delete w, psi, v;
    if (config.master) cout << "Memory deleted!\n";
}
