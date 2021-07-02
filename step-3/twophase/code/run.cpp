#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), fec(NULL), fespace(NULL),
    x_T(NULL), X_T(NULL),
    oper_T(NULL),
    ode_solver(NULL), cvode(NULL), arkode(NULL),
    paraview_out(NULL)
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
    delete m, k, t;
}

Flow_Operator::~Flow_Operator(){
    delete f, g;
    delete m, d, c, ct;
    delete M, D, C, Ct;
    delete A;
    delete psi, w;
}

Artic_sea::~Artic_sea(){
    delete pmesh, fec, fespace;
    delete x_T, X_T;
    delete oper_T;
    delete ode_solver;
    delete paraview_out;
    delete psi;
    if (config.master) cout << "Memory deleted \n";
}
