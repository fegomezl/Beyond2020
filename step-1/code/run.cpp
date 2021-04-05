#include "header.h"

Artic_sea::Artic_sea(bool master, int order, int refinements, bool last):
    master(master),
    order(order),
    refinements(refinements),
    last(last)
{}

void Artic_sea::run(const char *mesh_file){
    //Run the program
    make_grid(mesh_file);
    assemble_system();
    solve_system();
    output_results();
}

Artic_sea::~Artic_sea(){
    //Delete used memory
    delete pmesh;
    if (delete_fec) delete fec;
    delete fespace;
    delete a;
    delete b;
    delete x;
    delete u;
    if (master) cout << "Memory deleted!\n";
}
