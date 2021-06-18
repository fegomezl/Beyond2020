#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    p_exact_coeff(p_exact),
    v_exact_coeff(2, v_exact),
    pmesh(NULL),
    fec_rt(NULL), fec_l2(NULL),
    fespace_rt(NULL), fespace_l2(NULL),
    block_offsets(3), block_true_offsets(3),
    f(NULL), g(NULL),
    m(NULL), c(NULL),
    M(NULL), C(NULL), Ct(NULL),
    A(NULL),
    v(NULL), p(NULL)
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
    delete fec_rt, fec_l2;
    delete fespace_rt, fespace_l2;
    delete f, g;
    delete m, c;
    delete M, C, Ct;  
    delete A;
    delete v, p;
    if (config.master) cout << "Memory deleted!\n";
}
