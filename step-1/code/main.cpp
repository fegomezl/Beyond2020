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
    int order = 1;
    int refinements = 6;

    //Make program parameters readeable in execution
    OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&refinements, "-r", "--refinements",
                  "Number of total uniform refinements");

    //Check if parameters were read correctly
    args.Parse();
    if (!args.Good()){
      if (master) args.PrintUsage(cout);
      MPI_Finalize();
      return 1;
    }
    if (master) args.PrintOptions(cout);

    //Read mesh (serial)
    Mesh *mesh = new Mesh("data/star.mesh", 1, 1);
    int dim = mesh->Dimension();

    //Refine mesh (serial)
    int serial_refinements = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim);
    for (int ii = 0; ii < min(refinements, serial_refinements); ii++)
        mesh->UniformRefinement();

    //Make mesh (parallel), delete the serial
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    //Refine mesh (parallel)
    for (int ii = 0; ii < refinements - serial_refinements; ii++)
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
    } else {
        fec = new H1_FECollection(order = 1, dim);
        delete_fec = true;
    }
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_Int size = fespace->GlobalTrueVSize();

    //Set the boundary values to Dirichlet
    Array<int> ess_tdof_list;
    if (pmesh->bdr_attributes.Size()){
        Array<int> ess_bdr(pmesh->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    //Define biliniar form
    ParBilinearForm *a = new ParBilinearForm(fespace);
    ConstantCoefficient one(1.);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));
    a->Assemble();

    //Create RHS
    ParLinearForm *b = new ParLinearForm(fespace);
    b->AddDomainIntegrator(new DomainLFIntegrator(one));
    b->Assemble();

    //Define solution x
    ParGridFunction *x = new ParGridFunction(fespace);
    *x = 0.;

    //Create the linear system Ax=B
    HypreParMatrix A;
    Vector B, X;
    a->FormLinearSystem(ess_tdof_list, *x, *b, A, X, B);

    //Set the preconditioner
    HypreBoomerAMG *amg = new HypreBoomerAMG(A);
    amg->SetPrintLevel(0);

    //Solve the linear system Ax=B
    HyprePCG *pcg = new HyprePCG(A);
    pcg->SetPreconditioner(*amg);
    pcg->SetPrintLevel(0);
    pcg->SetTol(1e-12);
    pcg->SetMaxIter(200);
    pcg->Mult(B, X);

    //Recover the solution on each proccesor
    a->RecoverFEMSolution(X, *b, *x);

    //Print general information of the program
    if (master){ 
        cout << "Size: " << size << "\n"
             << "Serial refinements: " << serial_refinements << "\n"
             << "Parallel refinements: " << max(refinements-serial_refinements,0) << "\n"
             << "Total refinements: " << refinements << "\n";
    }

    //Output to Paraview
    ParaViewDataCollection paraview_out("graph", pmesh);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.SetCycle(0);
    paraview_out.SetTime(0.);
    paraview_out.SetLevelsOfDetail(order);
    paraview_out.RegisterField("Temperature", x);
    paraview_out.Save();

    //Delete used memory if needed
    delete pmesh;
    if (delete_fec) delete fec;
    delete fespace;
    delete a;
    delete b;
    delete x;
    delete amg;
    delete pcg;

    MPI_Finalize();
   
    return 0;
}
