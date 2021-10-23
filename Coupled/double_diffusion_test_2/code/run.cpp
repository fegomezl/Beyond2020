#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    dt(config.dt_init), last(false),
    vis_steps(config.vis_steps_max), vis_impressions(0),
    pmesh(NULL), fec(NULL), fespace(NULL),
    block_true_offsets(3),
    theta(NULL), phi(NULL), phase(NULL),
    cond_oper(NULL),
    ode_solver(NULL), arkode(NULL),
    paraview_out(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
        time_step();
    total_time = toc();
    output_results();
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete M;
    delete T;
    delete M_0;
    delete K_0;
    delete T_e;
    delete M_solver;
    delete T_solver;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec;
    delete fespace;
    delete theta;
    delete phi;
    delete phase;
    delete cond_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
