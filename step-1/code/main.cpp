#include <iostream>
#include <fstream>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]){

    //Define MPI parameters
    int nproc = 0, pid = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    bool master = (pid == 0);

    //Define program paramenters
    const char *mesh_file = "data/star.mesh";
    const char *device_config = "cpu";
    int order = 1;
    bool static_cond = false;
    bool pa = false;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", 
                   "-no-sc", "--no-static-condensation", 
                   "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", 
                   "-no-pa", "--no-partial-assembly", 
                   "Enable partial assembly.");
    args.AddOption(&device_config, "-d", "--device",
                   "Device configuration.");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
      if (master) args.PrintUsage(cout);
      MPI_Finalize();
      return 1;
    }
    if (master) args.PrintOptions(cout);

    //Set device (cpu, gpu...)
    Device device(device_config);
    if (master) device.Print();

    //Read mesh (serial)
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    //Refine mesh (serial)
    int ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
    for (int ii = 0; ii < ref_levels; ii++)
        mesh->UniformRefinement();

    //Make mesh (parallel), delete the serial
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    //Refine mesh (parallel)
    ref_levels = 2;
    for (int ii = 0; ii < ref_levels; ii++)
        pmesh->UniformRefinement();

    //Create the FEM space associated with the mesh
    FiniteElementCollection *fec;
    bool delete_fec;
    if (order > 0) {
        fec = new H1_FECollection(order, dim);
        delete_fec = true;
    } else if (pmesh->GetNodes()){
        fec = pmesh->GetNodes()->OwnFEC();
        delete_fec = false;
        if (master) 
            cout << "Using isoparamentric FEs: " << fec->Name() << "\n";
    } else {
        fec = new H1_FECollection(order = 1, dim);
        delete_fec = true;
    }
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_Int size = fespace->GlobalTrueVSize();
    if (master) cout << "Number of finite element unknowns: " << size << "\n";

    //Set the boundary values to Dirichlet
    Array<int> ess_tdof_list;
    if (pmesh->bdr_attributes.Size()){
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    //Create RHS
    ParLinearForm b(fespace);
    ConstantCoefficient one(1.);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble(); //Calculates values

    //Define solution x
    ParGridFunction x(fespace);
    x = 0.;

    //Define biliniar form
    ParBilinearForm a(fespace);
    if (pa) a.SetAssemblyLevel(AssemblyLevel::PARTIAL);  //Partial Assembly
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    if (static_cond) a.EnableStaticCondensation();       //Static Condensation
    a.Assemble();

    //Create the linear system
    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    //Solve the linear system AX=B
    Solver *prec = NULL;
    if (pa && UsesTensorBasis(*fespace))
        prec = new OperatorJacobiSmoother(a, ess_tdof_list);
    else 
        prec = new HypreBoomerAMG;
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    if (prec) cg.SetPreconditioner(*prec);
    cg.SetOperator(*A);
    cg.Mult(B, X);

    //Recover the solution on each proccesor
    a.RecoverFEMSolution(X, b, x);

    //Output to Paraview
    ParaViewDataCollection paraview_out("graph", pmesh);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.SetCycle(0);
    paraview_out.SetTime(0.0);
    paraview_out.SetLevelsOfDetail(order);
    paraview_out.RegisterField("Temperature", &x);
    paraview_out.Save();

    //Delete used memory if needed
    delete prec;
    delete pmesh;
    delete fespace;
    if (delete_fec) delete fec;

    MPI_Finalize();
   
    return 0;
}
