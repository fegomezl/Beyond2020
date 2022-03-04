#include "header.h"

Config::Config(int pid, int nproc):
    pid(pid),
    master(pid == 0),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    t(config.t_init), dt(config.dt_init), last(false),
    vis_steps(config.vis_steps_max), vis_impressions(0),
    pmesh(NULL), fec(NULL), fec_v(NULL), fespace(NULL), fespace_v(NULL),
    block_true_offsets(3),
    theta(NULL), phi(NULL), w(NULL), psi(NULL), v(NULL), rv(NULL), phase(NULL), 
    rV(NULL), V(NULL),
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
    delete M0_o;
    delete M1_o; 
    delete M0_e;  
    delete M1_e;  
    delete K0;
    delete K1;
    delete T0;
    delete T1;
    delete T0_e;
    delete T1_e;  
}

Flow_Operator::~Flow_Operator(){
    delete A00;
    delete A00_e;
    delete A10;
    delete A01_e;
    delete A01;
    delete A10_e;
    delete A11;
    delete A11_e;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec;
    delete fec_v;
    delete fespace;
    delete fespace_v;
    delete theta;
    delete phi;
    delete w;
    delete psi;
    delete v;
    delete rv;
    delete phase;
    delete rV;
    delete V;
    delete transport_oper;
    delete flow_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
