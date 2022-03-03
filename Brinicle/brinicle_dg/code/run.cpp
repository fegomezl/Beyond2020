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
    pmesh(NULL), fec(NULL), fec_dg(NULL), fec_v(NULL), fespace(NULL), fespace_dg(NULL), fespace_v(NULL),
    block_true_offsets(3),  block_true_offsets_dg(3),
    theta(NULL), phi(NULL), w(NULL), psi(NULL), v(NULL), rv(NULL), phase(NULL), 
    rV(NULL), V(NULL),
    cond_oper(NULL), flow_oper(NULL),
    ode_solver(NULL), arkode(NULL),
    paraview_out(NULL)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    if (config.t_final != 0)
        for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
            time_step();
    else
        iteration = vis_iteration = 0;
    total_time = toc();
    output_results();
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete M_theta; 
    delete M_phi; 
    delete K_theta; 
    delete K_phi; 
    delete T_theta; 
    delete T_phi; 
    delete B_theta;
    delete B_phi;
    delete B_dt_theta;
    delete B_dt_phi;
}

Flow_Operator::~Flow_Operator(){
    delete M;
    delete M_e;
    delete D;
    delete D_e;
    delete C;
    delete C_e;
    delete Ct;
    delete Ct_e;
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
    delete cond_oper;
    delete flow_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
