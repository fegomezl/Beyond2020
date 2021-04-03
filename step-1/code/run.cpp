#include "header.h"

Artic_sea::Artic_sea(bool master, int order, int refinements):
    master(master),
    order(order),
    refinements(refinements)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    solve_system();
    output_results();

    //Delete used memory
    delete pmesh;
    if (delete_fec) delete fec;
    delete fespace;
    delete a;
    delete b;
    delete x;
    delete u;
}
