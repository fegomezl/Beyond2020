// Gmsh project created on Tue Apr 20 21:20:33 2021
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;

Inner_Radius = 5;
Radius = 20;
Height = 20;

SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, -Height/2, 0, 0, Height, Radius, 2*Pi};
//+
Cylinder(2) = {0, 0, -Height/2, 0, 0, Height, Inner_Radius, 2*Pi};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Periodic Surface {7} = {6} Translate{0, 0, -Height};
//+
Physical Surface("Top", 1) = {6};
//+
Physical Surface("Bottom", 2) = {7};
//+
Physical Surface("Side", 3) = {5};
//+
Physical Surface("Inner_Side", 4) = {4};
//+
Physical Volume("Sea", 10) = {1};
