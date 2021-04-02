// Gmsh project created on Wed Mar 31 18:53:23 2021
SetFactory("OpenCASCADE");
//+
height = DefineNumber[ 20, Name "Parameters/height" ];
//+
Radius = DefineNumber[ 20, Name "Parameters/Radius" ];
//+
InnerRadius = DefineNumber[ 5, Name "Parameters/InnerRadius" ];
//+
Cylinder(1) = {0, 0, 0, 0, 0, height, Radius, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, height, InnerRadius, 2*Pi};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Surface{4} In Volume{1};
//+
Physical Surface("top", 10) = {6};
//+
Physical Surface("bottom", 11) = {7};
//+
Physical Surface("side", 12) = {5};
//+
Physical Surface("innerSide", 13) = {4};
//+
Physical Volume("Sea", 14) = {1};
