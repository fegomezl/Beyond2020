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
    block_true_offsets(3),
    w(NULL), psi(NULL), v(NULL),
    w_aux(NULL), psi_aux(NULL), v_aux(NULL), theta(NULL),
    g(NULL), f(NULL),
    m(NULL), d(NULL), c(NULL), 
    M(NULL), D(NULL), C(NULL)
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
    delete fec_v;
    delete fespace;
    delete fespace_v;
    delete w;
    delete psi;
    delete v;
    delete w_aux;
    delete psi_aux;
    delete v_aux;
    delete theta;
    delete g;
    delete f;
    delete m;
    delete d;
    delete c;
    delete M;
    delete D;
    if (config.master) cout << "Memory deleted!\n";
}
