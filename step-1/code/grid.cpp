#include "header.h"

void Artic_sea::make_grid(){
    //Read mesh (serial)
    Mesh *mesh = new Mesh("data/try2.msh", 1, 1);
    dim = mesh->Dimension();

    //Refine mesh (serial)
    serial_refinements = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
    for (int ii = 0; ii < min(refinements, serial_refinements); ii++)
        mesh->UniformRefinement();

    //Make mesh (parallel), delete the serial
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    //Refine mesh (parallel)
    for (int ii = 0; ii < refinements - serial_refinements; ii++)
        pmesh->UniformRefinement();

    //Create the FEM space associated with the mesh
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
        delete_fec = true;
    } else if (pmesh->GetNodes()){
        fec = pmesh->GetNodes()->OwnFEC();
        delete_fec = false;
    } else {
        fec = new H1_FECollection(order = 1, dim);
        delete_fec = true;
    }
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();
}
