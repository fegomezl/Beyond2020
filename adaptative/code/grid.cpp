#include "header.h"

void Artic_sea::make_grid(const char *mesh_file){
    //Read mesh (serial)
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    dim = mesh->Dimension();
    int sdim = mesh->SpaceDimension();//**************************
    
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
    
    //*********************************
       MFEM_VERIFY(pmesh.bdr_attributes.Size() > 0,
               "Boundary attributes required in the mesh.");
   Array<int> ess_bdr(pmesh.bdr_attributes.Max());
   ess_bdr = 1;
   //**********************************

    //Calculate minimum size of elements
    double null;
    pmesh->GetCharacteristics(h_min, null, null, null);

    //Create the FEM space associated with the mesh
    if (config.order > 0)
        fec = new H1_FECollection(config.order, dim);
    else
        fec = new H1_FECollection(config.order = 1, dim);
    fespace = new ParFiniteElementSpace(pmesh, fec);
    size = fespace->GlobalTrueVSize();

   //**********************************************
   ParBilinearForm a(&fespace);
   ParLinearForm b(&fespace);

   ConstantCoefficient one(1.0);
   FunctionCoefficient bdr(bdr_func);
   FunctionCoefficient rhs(rhs_func);

   BilinearFormIntegrator *integ = new DiffusionIntegrator(one);
   a.AddDomainIntegrator(integ);
   b.AddDomainIntegrator(new DomainLFIntegrator(rhs));

   // 8. The solution vector x and the associated finite element grid function
   //    will be maintained over the AMR iterations.
   ParGridFunction x(&fespace);
}
