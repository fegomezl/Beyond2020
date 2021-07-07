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

    //Refine mesh (parallel)
    for (int ii = 0; ii < config.refinements - serial_refinements; ii++)
        pmesh->UniformRefinement();

    //Calculate minimum size of elements
    double null;
    pmesh->GetCharacteristics(h_min, null, null, null);

    //Create the FEM spaces associated with the mesh
    fec = new H1_FECollection(config.order, dim);
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();

    fec_v = new RT_FECollection(config.order, dim);
    fespace_v = new ParFiniteElementSpace(pmesh, fec_v);
    size_v = fespace_v->GlobalTrueVSize();
    
    //Create the block offsets
    block_offsets[0] = 0;
    block_offsets[1] = fespace->GetVSize();
    block_offsets[2] = fespace->GetVSize();
    block_offsets.PartialSum();

    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace->TrueVSize();
    block_true_offsets[2] = fespace->TrueVSize();
    block_true_offsets.PartialSum();

    //Initialize the corresponding vectors
    x.Update(block_offsets); X.Update(block_true_offsets);
    b.Update(block_offsets); B.Update(block_true_offsets);
}
