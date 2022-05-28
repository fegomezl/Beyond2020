#include "header.h"

//Create the mesh and the FES
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

    if (config.master)
        cout << "\nSize (H1): " << size_H1 << "\n"
             << "Size (ND): " << size_ND << "\n"
             << "Mesh Size: " << h_min << " (" << h_min*L_ref << " mm)\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << config.refinements - serial_refinements << "\n"
             << "Total refinements: " << config.refinements << "\n\n";

    //Create the FEM spaces associated with the mesh
    fec_H1 = new H1_FECollection(config.order, dim);
    fespace_H1 = new ParFiniteElementSpace(pmesh, fec_H1);
    size_H1 = fespace_H1->GlobalTrueVSize();

    fec_ND = new ND_FECollection(config.order, dim);
    fespace_ND = new ParFiniteElementSpace(pmesh, fec_ND);
    size_ND = fespace_ND->GlobalTrueVSize();

    //Create the block offsets
    block_offsets_H1[0] = 0;
    block_offsets_H1[1] = fespace_H1->TrueVSize();
    block_offsets_H1[2] = fespace_H1->TrueVSize();
    block_offsets_H1.PartialSum();

    X.Update(block_offsets_H1); X = 0.;
    Y.Update(block_offsets_H1); Y = 0.;
}
