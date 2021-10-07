#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), fec(NULL), fespace(NULL),
    x(NULL), X(NULL),
    oper(NULL),
    ode_solver(NULL), cvode(NULL), arkode(NULL),
    exact(initial),
    paraview_out(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    //for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
    //    time_step();
    total_time = toc();
    output_results();
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete m;
    delete k;
    delete M;
    delete M_e;
    delete M_0;
    delete K_0;
    delete T;
    delete T_e;
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
    if (config.master) cout << "Memory deleted \n";
}
