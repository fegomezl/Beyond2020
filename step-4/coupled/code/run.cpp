#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), fec(NULL), fec_v(NULL), fespace(NULL), fespace_v(NULL),
    block_true_offsets(3),
    theta(NULL), phi(NULL), phase(NULL), Theta(NULL),
    oper_T(NULL),
    ode_solver(NULL), cvode(NULL), arkode(NULL),
    paraview_out(NULL), flow_oper(NULL), x_psi(NULL), x_v(NULL), x_w(NULL)
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
    delete m_theta;
    delete m_phi;
    delete k_theta;
    delete k_phi;
    delete t_theta;
    delete t_phi;
}

Flow_Operator::~Flow_Operator(){
    delete f;
    delete g;
    delete m;
    delete d;
    delete c;
    delete M;
    delete D;
    delete psi;
    delete w;
    delete v;
    delete psi_aux;
    delete w_aux;
    delete v_aux;
    delete theta;
   // delete H;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec;
    delete fespace;
    delete fec_v;
    delete fespace_v;
    delete theta;
    delete phi;
    delete phase;
    delete Theta;
    delete oper_T;
    delete ode_solver;
    delete paraview_out;
    delete flow_oper;
    delete x_psi;
    delete x_v;
    delete x_w;
    if (config.master) cout << "Memory deleted \n";
}
