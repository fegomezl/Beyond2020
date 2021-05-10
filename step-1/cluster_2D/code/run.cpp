#include "header.h"

Config::Config(bool master, int nproc):
    master(master), 
    nproc(nproc)
{}

Artic_sea::Artic_sea(Config config):
    config(config),
    u(exact),
    r(rf),
    pmesh(NULL),
    fec(NULL),
    fespace(NULL),
    a(NULL),
    b(NULL),
    x(NULL)
{}

void Artic_sea::run(const char *mesh_file, double &H_min, double &Total_time, double &Mean_error){
    //Run the program
    tic();
    make_grid(mesh_file);
    assemble_system();
    solve_system();
    Total_time += toc();
    H_min += h_min;
    Mean_error += l2_error;
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    delete fec;
    delete fespace;
    delete a;
    delete b;
    delete x;
}
