#include "header.h"

Config::Config(int pid, int nproc):
    pid(pid),
    master(pid == 0),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    t(config.t_init), dt(config.dt_init), last(false),
    vis_steps(config.vis_steps_max), vis_print(0),
    pmesh(NULL), 
    fec_H1(NULL), fec_L2(NULL), fec_ND(NULL), 
    fespace_H1(NULL), fespace_L2(NULL), fespace_ND(NULL),
    block_offsets_H1(3),
    block_offsets_L2(3),
    temperature(NULL), salinity(NULL), phase(NULL), 
    vorticity(NULL), stream(NULL), 
    velocity(NULL), rvelocity(NULL), 
    Velocity(NULL), rVelocity(NULL),
    transport_oper(NULL), flow_oper(NULL),
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

Transport_Operator::~Transport_Operator(){
    //Delete used memory
    delete M0; 
    delete M1; 
    delete K0; 
    delete K1; 
    delete T0; 
    delete T1; 
    delete B0;
    delete B1;
}

Flow_Operator::~Flow_Operator(){
    delete B0;
    delete B1;
    delete A00;
    delete A01;
    delete A10;
    delete A11;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec_H1;
    delete fec_L2;
    delete fec_ND;
    delete fespace_H1;
    delete fespace_L2;
    delete fespace_ND;
    delete temperature;
    delete salinity;
    delete phase;
    delete vorticity;
    delete stream;
    delete velocity;
    delete rvelocity;
    delete Velocity;
    delete rVelocity;
    delete transport_oper;
    delete flow_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
