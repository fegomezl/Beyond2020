// Gmsh project created on Thu Sep 16 23:50:02 2021
SetFactory("OpenCASCADE");
Rmin = 0;
Rmax = 10;
Zmin = 0;
Zmax = 10;
Hside = 2;

NR = 30;
NZ = 30;

Point(1) = {Rmin, Zmin, 0, 1.0};
//+
Point(2) = {Rmax, Zmin, 0, 1.0};
//+
Point(3) = {Rmax, Zmax, 0, 1.0};
//+
Point(4) = {Rmin, Zmax, 0, 1.0};
//+
Point(5) = {Rmin, (Zmax-Zmin+Hside)/2, 0, 1.0};
//+
Point(6) = {Rmin, (Zmax-Zmin-Hside)/2, 0, 1.0};
//+
Point(7) = {Rmax, (Zmax-Zmin+Hside)/2, 0, 1.0};
//+
Point(8) = {Rmax, (Zmax-Zmin-Hside)/2, 0, 1.0};
//+
Point(9) = {(Rmax-Rmin-Hside)/2, (Zmax-Zmin+Hside)/2, 0, 1.0};
//+
Point(10) = {(Rmax-Rmin-Hside)/2, (Zmax-Zmin-Hside)/2, 0, 1.0};
//+
Point(11) = {(Rmax-Rmin+Hside)/2, (Zmax-Zmin+Hside)/2, 0, 1.0};
//+
Point(12) = {(Rmax-Rmin+Hside)/2, (Zmax-Zmin-Hside)/2, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {3, 7};
//+
Line(5) = {6, 1};
//+
Line(6) = {8, 2};
//+
Line(7) = {5, 9};
//+
Line(8) = {11, 7};
//+
Line(9) = {6, 10};
//+
Line(10) = {12, 8};
//+
Line(11) = {9, 10};
//+
Line(12) = {11, 12};
//+
//+
Curve Loop(1) = {2, 3, 7, 11, -9, 5, 1, -6, -10, -12, 8, -4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Down", 1) = {1};
//+
Physical Curve("Top", 2) = {2};
//+
Physical Curve("Left", 3) = {3, 5};
//+
Physical Curve("Hleft", 4) = {7, 11, 9};
//+
Physical Curve("Right", 5) = {4, 6};
//+
Physical Curve("Hright", 6) = {8, 12, 10};
//+
Physical Surface("Sea", 7) = {1};
