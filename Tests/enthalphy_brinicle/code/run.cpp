#include "header.h"

Config::Config(bool master, int nproc):
    master(master),
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    pmesh(NULL), 
    fec_H1(NULL),     fec_ND(NULL),
    fespace_H1(NULL), fespace_ND(NULL),
    block_offsets_H1(3),
    enthalphy(NULL), salinity(NULL), 
    vorticity(NULL), stream(NULL), 
    velocity(NULL), r_velocity(NULL), 
    Velocity(NULL), r_Velocity(NULL),
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
    delete m0;
    delete m1;
    delete k0;
    delete k1;
    delete b0;
    delete b1;
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
    delete T0_e;
    delete B0;
    delete B1;
    delete B0_dt;
    delete B1_dt;
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
    delete fec_H1;
    delete fec_ND;
    delete fespace_H1;
    delete fespace_ND;
    delete enthalphy;
    delete salinity;
    delete vorticity;
    delete stream;
    delete velocity;
    delete r_velocity;
    delete Velocity;
    delete r_Velocity;
    delete transport_oper;
    delete flow_oper;
    delete ode_solver;
    delete paraview_out;
    if (config.master) cout << "Memory deleted \n";
}
