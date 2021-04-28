//+
SetFactory("OpenCASCADE");
R = 10;
Z = 10;

NR = 10;
NZ = 10;

Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {R, 0, 0, 1.0};
//+
Point(3) = {R, Z, 0, 1.0};
//+
Point(4) = {0, Z, 0, 1.0};
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
Transfinite Curve {3, 1} = NR Using Progression 1;
//+
Transfinite Curve {4, 2} = NZ Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
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
