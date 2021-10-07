//+
SetFactory("OpenCASCADE");
Rmin = 0;
Rmax = 25;
Zmin = 0;
Zmax = 25;

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
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("down", 1) = {1};
//+
Physical Curve("top", 2) = {3};
//+
Physical Curve("left", 3) = {4};
//+
Physical Curve("right", 4) = {2};
//+
Physical Surface("Sea", 5) = {1};
