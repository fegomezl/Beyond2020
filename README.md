# Project Beyond 2020-I

## Running the programs

1. Copy the hidden file **.template\_local\_config.txt** with the name **local\_config.mk**
2. Add the MFEM instalation directory.
3. If you want to move the graphs to other folder after running a program, change the **NULL** option of the **SHARE\_DIR** variable to the folder directory.

## List of programs

- Step 1
  - 3D:
    Solution of a benchmark Poisson equation problem with convergence analysis on a 3D hollow cylinder.
  - 2D:
    Solution of a benchmark Poisson equation problem with convergence analysis on a 2D square which simulates a hollow cylinder on cylindrical coordinates.
- Step 2
  - onephase:
    Solution of a benchmark time dependent heat equation with convergence analysis on a 2D square which simulates a hollow cylinder on cylindrical coordinates.
  - twophase:
    Simulation of a Stefan problem in a 2D square which simulates a cylinder on cylindrical coordinates.
