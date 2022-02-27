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
    fec_H1 = new H1_FECollection(config.order, dim);
    fec_L2 = new L2_FECollection(config.order, dim);
    fec_ND = new ND_FECollection(config.order, dim);

    fespace_H1 = new ParFiniteElementSpace(pmesh, fec_H1);
    fespace_L2 = new ParFiniteElementSpace(pmesh, fec_L2);
    fespace_ND = new ParFiniteElementSpace(pmesh, fec_ND);

    size_H1 = fespace_H1->GlobalTrueVSize();
    size_L2 = fespace_L2->GlobalTrueVSize();
    size_ND = fespace_ND->GlobalTrueVSize();

    //Create the block offsets
    block_offsets_H1[0] = 0;
    block_offsets_H1[1] = fespace_H1->TrueVSize();
    block_offsets_H1[2] = fespace_H1->TrueVSize();
    block_offsets_H1.PartialSum();

    block_offsets_L2[0] = 0;
    block_offsets_L2[1] = fespace_L2->TrueVSize();
    block_offsets_L2[2] = fespace_L2->TrueVSize();
    block_offsets_L2.PartialSum();

    //Update block vectors
    X.Update(block_offsets_L2);
    Y.Update(block_offsets_H1);
}
