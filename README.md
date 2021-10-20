# Heat and Mass Transfer with Phase Change

In this repository you will find a list of simulation regarding the effect of phase change on heat and mass transfer. All of the programs use the Finite Element Method, which is implmented on C++ with the library MFEM.

In each of the simulations, the domain is a cylinder, sometimes hollow, with an axial symmetry. Inside the domain, it is assumed that the material is water, which can be liquid and solid. 

## Running the programs

1. Copy the hidden file **.template\_local\_config.txt** with the name **local\_config.mk**
2. Add the MFEM instalation directory.
3. If you want to move the graphs to other folder after running a program, change the **NULL** option of the **SHARE\_DIR** variable to the folder directory.
4. If you want to chage some quantity on a simulation, in each folder the file **settings/parameters.txt** has all the main parameters of the simulation.
5. To run a simulation, you only have to write **make**.

## List of programs

- Time_Independent
  - 3D:
    Benchmark Poisson equation problem with convergence analysis on a 3D hollow cylinder.
  - 2D:
    Benchmark Poisson equation problem with convergence analysis on a 2D square which simulates a hollow cylinder on cylindrical coordinates (the rest of the simulations assume this domain representation).
- Time_Dependent
  - onephase:
    Benchmark time dependent heat equation with convergence analysis, in which only one phase is present.
  - twophase:
    Simulation with initial conditions set on different phases.
  - twophase_exact:
    Stefan problem simulation with comparison according to the exact Stefan solution.
- Flow
  - vorticity:
    Quasi-static flow problem with a stream function-vorticity formulation and a frozen region.
  - vorticity_barrier:
    Quasi-static flow problem with a stream function-vorticity formulation and an obstacle.
  - vorticity_exact:
    Benchmark quasi-static flow problem with a stream function-vorticity formulation, including a convergence analysis.
- Coupled
  - heat_flow:
    Heat diffusion problem with convection effects, where the flow is quasi-static and approched by a stream function-vorticity formulation.
  - double_diffusion:
    Heat and salinity diffusion problem, where both variables are weakly coplued.
  - diffusion_flow:
    Heat and salinity diffusion problem, where both variables are weakly coplued, including convection effects with a quasi-static flow and a stream function-vorticity approach.
- Brinicle
  - Simulation from **diffusion_flow** applied to form a brinicle, which is an ice channel formed on the artic sea because of icy brine fluxes that enter the ocean from the top layer of ice.
