#include "header.h"

void Artic_sea::make_grid(const char *mesh_file){
    //Read mesh (serial)
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    dim = mesh->Dimension();

    //Refine mesh (serial)
    config.serial_refinements = min(config.refinements, config.serial_refinements);
    for (int ii = 0; ii < config.serial_refinements; ii++)
        mesh->UniformRefinement();

    //Make mesh (parallel), delete the serial
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    //Refine mesh (parallel)
    for (int ii = 0; ii < config.refinements - config.serial_refinements; ii++)
        pmesh->UniformRefinement();

    //Calculate minimum size of elements
    double null;
    pmesh->GetCharacteristics(h_min, null, null, null);

    //Create the FEM space associated with the mesh
    if (config.order > 0) {
        fec = new H1_FECollection(config.order, dim);
        delete_fec = true;
    } else if (pmesh->GetNodes()){
        fec = pmesh->GetNodes()->OwnFEC();
        delete_fec = false;
    } else {
        fec = new H1_FECollection(config.order = 1, dim);
        delete_fec = true;
    }
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();
}
