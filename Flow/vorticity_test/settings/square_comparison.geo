//+
SetFactory("OpenCASCADE");
Rmin = 0;
Rmax = 10;
Zmin = 0;
Zmax = 10;
Hside = 2;

NR = 100;
NZ = 100;

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
Line(13) = {5,6};
//+
Line(14) = {7,8};
//+
Curve Loop(1) = {2, 3, 7, 11, -9, 5, 1, -6, -10, -12, 8, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 11, -9, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 14, -10, -12};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {7, 9, 8, 10} = 5 Using Progression 1;
//+
Transfinite Curve {13, 11, 12, 14} = 3 Using Progression 1;
//+
Transfinite Surface {2} = {5, 9, 10, 6};
//+
Transfinite Surface {3} = {11, 7, 8, 12};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Physical Curve("top", 1) = {2};
//+
Physical Curve("bottom", 2) = {1};
//+
Physical Curve("left", 3) = {3, 13, 5};
//+
Physical Curve("right", 4) = {4, 14, 6};
//+
Physical Surface("Sea", 5) = {2, 1, 3};

