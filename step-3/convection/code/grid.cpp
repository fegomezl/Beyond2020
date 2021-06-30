#include "header.h"

void Artic_sea::make_grid(const char *mesh_file){
    //Read mesh (serial)
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    dim = mesh->Dimension();

    //Calculate how many serial refinements are needed
    //More than 1000 cells per processor
    int elements = mesh->GetNE();
    int min_elements = 1000.*config.nproc;
    if (min_elements > elements)
        serial_refinements = min(config.refinements, (int)floor(log(min_elements/elements)/(dim*log(2.))));
    else
        serial_refinements = 0;

    //Refine mesh (serial)
    for (int ii = 0; ii < serial_refinements; ii++)
        mesh->UniformRefinement();

    //Make mesh (parallel), delete the serial
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    //Calculate minimum size of elements
    double null;
    pmesh->GetCharacteristics(h_min, null, null, null);

    //Refine mesh (parallel)
    for (int ii = 0; ii < config.refinements - serial_refinements; ii++)
        pmesh->UniformRefinement();

    //Create the FEM space associated with the mesh
    if (config.order > 0)
        fec = new H1_FECollection(config.order, dim);
    else
        fec = new H1_FECollection(config.order = 1, dim);
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();

    //Create the FEM space associated with the mesh
    fec_psi = new H1_FECollection(config.order, dim);
    fespace_psi = new ParFiniteElementSpace(pmesh, fec_psi);

    fec_w = new H1_FECollection(config.order, dim);
    fespace_w = new ParFiniteElementSpace(pmesh, fec_w);

    fec_v = new RT_FECollection(config.order, dim);
    fespace_v = new ParFiniteElementSpace(pmesh, fec_v);
}
