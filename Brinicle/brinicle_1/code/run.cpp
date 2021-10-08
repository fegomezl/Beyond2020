#include "header.h"

Config::Config(int pid, int nproc):
    pid(pid),
    master(pid == 0),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), fec(NULL), fec_v(NULL), fespace(NULL), fespace_v(NULL),
    block_true_offsets(3),
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
    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
        time_step();
    total_time = toc();
    output_results();
}

Conduction_Operator::~Conduction_Operator(){
    //Delete used memory
    delete m_theta; 
    delete m_phi;
    delete k_theta;
    delete k_phi;
    delete f_theta;
    delete f_phi;
    delete M_theta; 
    delete M_e_theta; 
    delete M_0_theta; 
    delete M_phi; 
    delete M_e_phi; 
    delete M_0_phi; 
    delete K_0_theta; 
    delete K_0_phi; 
    delete T_theta; 
    delete T_e_theta; 
    delete T_phi; 
    delete T_e_phi; 
}

Flow_Operator::~Flow_Operator(){
    delete f;
    delete g;
    delete m;
    delete d;
    delete c;
    delete ct;
    delete M;
    delete D;
    delete C;
    delete Ct;
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
