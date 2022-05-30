SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthFactor=r;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {R, 0, 0, 1.0};
Point(3) = {R, Z, 0, 1.0};
Point(4) = {0, Z, 0, 1.0};
Point(5) = {r, Z, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 4};
Line(5) = {4, 1};
Curve Loop(1) = {5, 1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Curve("Down", 1)  = {1};
Physical Curve("Top", 2)   = {3};
Physical Curve("Left", 3)  = {5};
Physical Curve("Right", 4) = {2};
Physical Curve("Inlet", 5) = {4};
Physical Surface("Sea", 6) = {1};
