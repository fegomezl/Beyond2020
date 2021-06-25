#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include "stokes.hpp"

// example runs
// mpirun -np 4 ./stokes -m ./ball2D.msh -petscopts ./stokes_fieldsplit
// mpirun -np 4 ./stokes -m ./ball2D.msh -petscopts ./stokes_fieldsplit_01

int main(int argc, char *argv[])
{
   // Initialize MPI.
   int nprocs, myrank;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);



   // Parse command-line options.
   const char *mesh_file = "../../data/star.mesh";
   int order = 1;
   bool static_cond = false;
   int ser_ref_levels = 1;
   int par_ref_levels = 1;
   double newton_rel_tol = 1e-7;
   double newton_abs_tol = 1e-12;
   int newton_iter = 10;
   int print_level = 1;
   bool visualization = false;

   const char *petscrc_file = "stokes_fieldsplit";

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ser_ref_levels,
                  "-rs",
                  "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&par_ref_levels,
                  "-rp",
                  "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&visualization,
                  "-vis",
                  "--visualization",
                  "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&newton_rel_tol,
                  "-rel",
                  "--relative-tolerance",
                  "Relative tolerance for the Newton solve.");
   args.AddOption(&newton_abs_tol,
                  "-abs",
                  "--absolute-tolerance",
                  "Absolute tolerance for the Newton solve.");
   args.AddOption(&newton_iter,
                  "-it",
                  "--newton-iterations",
                  "Maximum iterations for the Newton solve.");
   args.AddOption(&petscrc_file, "-petscopts", "--petscopts",
                                 "PetscOptions file to use.");
   args.Parse();
   if (!args.Good())
   {
      if (myrank == 0)
      {
         args.PrintUsage(std::cout);
      }
      MPI_Finalize();
      return 1;
   }

   if (myrank == 0)
   {
      args.PrintOptions(std::cout);
   }

   mfem::MFEMInitializePetsc(NULL,NULL,petscrc_file,NULL);

   // Read the (serial) mesh from the given mesh file on all processors.  We
   // can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   // and volume meshes with the same code.
   mfem::Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // Refine the serial mesh on all processors to increase the resolution. In
   // this example we do 'ref_levels' of uniform refinement. We choose
   // 'ref_levels' to be the largest number that gives a final mesh with no
   // more than 10,000 elements.
   {
      int ref_levels =
         (int)floor(log(100./mesh.GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
   }

   // Define a parallel mesh by a partitioning of the serial mesh. Refine
   // this mesh further in parallel to increase the resolution. Once the
   // parallel mesh is defined, the serial mesh can be deleted.
   mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
   /*
   {
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh.UniformRefinement();
      }
   }*/


   mfem::StokesSolver* solver=new mfem::StokesSolver(&pmesh,2);

   mfem::ConstantCoefficient viscosity(1);
   solver->SetViscosity(viscosity);
   solver->AddVelocityBC(2,0,1.0);
   solver->AddVelocityBC(1,4,0.0);

   mfem::Vector vload(2); vload(0)=1.0; vload(1)=0.0; //vload(2)=0.0;
   mfem::VectorConstantCoefficient load(vload);
   solver->SetVolForces(load);

   //solver->AddVelocityBC(3,4,0.0);
   //solver->AddVelocityBC(1,4,0.0);

   //mfem::ConstantCoefficient bcvelx(2.0);
   //mfem::ConstantCoefficient bcvely(1.5);
   //solver->AddVelocityBC(1,0,bcvelx);
   //solver->AddVelocityBC(2,1,bcvely);

   solver->Solve();


   {
       mfem::ParGridFunction& veloc=solver->GetVelocity();
       mfem::ParGridFunction& press=solver->GetPressure();

       mfem::ParaViewDataCollection paraview_dc("Stokes", &pmesh);
       paraview_dc.SetPrefixPath("ParaView");
       paraview_dc.SetLevelsOfDetail(order);
       paraview_dc.SetDataFormat(mfem::VTKFormat::BINARY);
       paraview_dc.SetHighOrderOutput(true);
       paraview_dc.SetCycle(0);
       paraview_dc.SetTime(0.0);
       paraview_dc.RegisterField("velocity",&veloc);
       paraview_dc.RegisterField("pressure",&press);
       paraview_dc.Save();
   }



   delete solver;

   mfem::MFEMFinalizePetsc();
   MPI_Finalize();
   return 0;
}


