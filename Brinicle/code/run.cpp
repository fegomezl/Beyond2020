#include "header.h"

//Initialization of the program variables 
Config::Config(int pid, int nproc):
    pid(pid),
    master(pid == 0),
    nproc(nproc)
{}

//Initialization of the program
Artic_sea::Artic_sea(Config config):
    config(config),
    t(0.), dt(config.dt_init), last(false),
    vis_steps(config.vis_steps_max), vis_print(0),
    pmesh(NULL), 
    fec_H1(NULL), fespace_H1(NULL), 
    block_offsets_H1(3),
    transport_oper(NULL), flow_oper(NULL),
    ode_solver(NULL), arkode(NULL),
    paraview_out(NULL)
{}

//Run the program
void Artic_sea::run(const char *mesh_file){
    make_grid(mesh_file);
    assemble_system();
    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++)
        time_step();
    total_time = toc();

    if (config.master){
        cout << "\n\nIterations: " << iteration-1 << "\n"
             << "Total Printed: " << vis_print << "\n"
             << "Total execution time: " << total_time << " s" << "\n\n";
    }
    if (config.master){
        ofstream out;
        out.open("results/graph/settings.txt", std::ios::app);
        out << "\n\nIterations: " << iteration-1 << "\n"
             << "Total Printed: " << vis_print << "\n"
             << "Total execution time: " << total_time << " s" << "\n\n";
        out.close();
    }
}

//Delete used memory

Transport_Operator::~Transport_Operator(){
    delete M0; 
    delete M1; 
    delete K0; 
    delete K1; 
    delete T0; 
    delete T1; 
}

Flow_Operator::~Flow_Operator(){
    delete A00;
    delete A01;
    delete A10;
    delete A11;
}

Artic_sea::~Artic_sea(){
    delete pmesh;
    delete fec_H1;
    delete fespace_H1;
    delete transport_oper;
    delete flow_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
