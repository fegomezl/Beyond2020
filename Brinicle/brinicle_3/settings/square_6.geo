//+
SetFactory("OpenCASCADE");
Rmin = 0;
Rmax = 20;
Zmin = 0;
Zmax = 20;
l = 2;
h = 5;

NR = 10;
NZ = 10;

Point(1) = {Rmin, Zmin, 0, 1.0};
//+
Point(2) = {Rmax, Zmin, 0, 1.0};
//+
Point(3) = {Rmax, Zmax, 0, 1.0};
//+
Point(4) = {Rmin, Zmax, 0, 1.0};
//+
Point(5) = {l, Zmax, 0, 1.0};
//+
Point(6) = {Rmax, h, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 6};
//+
Line(3) = {6, 3};
//+
Line(4) = {3, 5};
//+
Line(5) = {5, 4};
//+
Line(6) = {4, 1};//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("down", 1) = {1};
//+
Physical Curve("top", 2) = {4};
//+
Physical Curve("left", 3) = {6};
//+
Physical Curve("right", 4) = {3};
//+
Physical Curve("l", 5) = {5};
//+
Physical Curve("h", 6) = {2};
//+
Physical Surface("Sea", 7) = {1};
