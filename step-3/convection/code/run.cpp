#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL),
    fec(NULL),
    fespace(NULL),
    initial_f(initial),
    x(NULL),
    X(NULL),
    oper(NULL),
    ode_solver(NULL),
    cvode(NULL),
    arkode(NULL),
    paraview_out(NULL),
    fec_psi(NULL), fec_w(NULL), fec_v(NULL),
    fespace_psi(NULL), fespace_w(NULL), fespace_v(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
        time_step();
    output_results();
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete m;
    delete k;
    delete t;
}

Flow_Operator::~Flow_Operator(){
    delete f, g;
    delete m, c, ct;
    delete M, C, Ct;
    delete A;
    delete w, psi, v;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec;
    delete fespace;
    delete x;
    delete X;
    delete oper;
    delete ode_solver;
    delete paraview_out;
    delete fec_psi, fec_w, fec_v;
    delete fespace_psi, fespace_w, fespace_v;
    if (config.master) cout << "Memory deleted \n";
}
