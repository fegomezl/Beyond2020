// Gmsh project created on Sun Mar 13 16:06:19 2022
SetFactory("OpenCASCADE");
//+
R=10;
Z=20;
Cylinder(1) = {0, 0, 0, 0, 0, Z, R, 2*Pi};
