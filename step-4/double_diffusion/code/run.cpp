#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), fec(NULL), fespace(NULL),
    block_true_offsets(3),
    theta(NULL), phi(NULL), phase(NULL),
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
    delete m_theta, m_phi;
    delete k_theta, k_phi;
    delete t_theta, t_phi;
}

Artic_sea::~Artic_sea(){
    delete pmesh, fec, fespace;
    delete theta, phi, phase;
    delete oper_T;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
