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

    //Create the FEM space associated with the mesh
    fec_rt = new RT_FECollection(config.order, dim);
    fespace_rt = new ParFiniteElementSpace(pmesh, fec_rt);
    size_rt = fespace_rt->GlobalTrueVSize();

    fec_l2 = new L2_FECollection(config.order, dim);
    fespace_l2 = new ParFiniteElementSpace(pmesh, fec_l2);
    size_l2 = fespace_l2->GlobalTrueVSize();
    
    //Create the block offsets
    block_offsets[0] = 0;
    block_offsets[1] = fespace_rt->GetVSize();
    block_offsets[2] = fespace_l2->GetVSize();
    block_offsets.PartialSum();

    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace_rt->TrueVSize();
    block_true_offsets[2] = fespace_l2->TrueVSize();
    block_true_offsets.PartialSum();

    //Initialize the corresponding vectors
    const char *device_config = "cpu";
    Device device(device_config);
    MemoryType mt = device.GetMemoryType();
    x.Update(block_offsets, mt); X.Update(block_true_offsets, mt);
    b.Update(block_offsets, mt); B.Update(block_true_offsets, mt);
}
