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
    
    if (!config.restart){
        //Refine mesh (serial)
        for (int ii = 0; ii < serial_refinements; ii++)
            mesh->UniformRefinement();
    }
    
    //Make mesh (parallel), delete the serial
    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    
    if (!config.restart){
        //Refine mesh (parallel)
        for (int ii = 0; ii < config.refinements - serial_refinements; ii++)
            pmesh->UniformRefinement();
    }

    if (config.restart){
        std::ifstream in;
        std::ostringstream oss;
        oss << std::setw(10) << std::setfill('0') << config.pid;
        std::string n_mesh = "results/restart/pmesh_"+oss.str()+".msh";

        in.open(n_mesh.c_str(),std::ios::in);
        pmesh->Load(in,1,0); 
        in.close();
    }

    //Calculate minimum size of elements
    double null;
    pmesh->GetCharacteristics(h_min, null, null, null);

    //Create the FEM space associated with the mesh
    fec = new H1_FECollection(config.order, dim);
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();

    fec_v = new ND_FECollection(config.order, dim);
    fespace_v = new ParFiniteElementSpace(pmesh, fec_v);
    size_v = fespace_v->GlobalTrueVSize();

    //Create the block offsets
    block_true_offsets[0] = 0;
    block_true_offsets[1] = fespace->TrueVSize();
    block_true_offsets[2] = fespace->TrueVSize();
    block_true_offsets.PartialSum();

    X.Update(block_true_offsets);
    Z.Update(block_true_offsets);
}
