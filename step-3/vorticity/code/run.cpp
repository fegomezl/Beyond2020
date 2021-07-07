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
    block_offsets(3), block_true_offsets(3),
    w(NULL), psi(NULL), v(NULL),
    g(NULL), f(NULL),
    m(NULL), d(NULL), c(NULL), 
    M(NULL), D(NULL), C(NULL),
    w_aux(NULL), psi_aux(NULL)
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
    delete fec, fec_v;
    delete fespace, fespace_v;
    delete w, psi, v;
    delete g, f;
    delete m, d, c;
    delete M, D, C;  
    delete w_aux, psi_aux;
    if (config.master) cout << "Memory deleted!\n";
}
