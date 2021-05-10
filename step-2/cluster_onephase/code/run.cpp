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
    arkode(NULL)
{}

void Artic_sea::run(const char *mesh_file, double &H_size, double &Total_time, double &Mean_error){
    tic();

    //Run the program
    make_grid(mesh_file);
    assemble_system();
    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
        time_step();

    Total_time += toc();
    H_size += h_size;
    Mean_error += total_error/vis_impressions;
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete m;
    delete k;
    delete T;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec;
    delete fespace;
    delete x;
    delete X;
    delete oper;
    delete ode_solver;
}
