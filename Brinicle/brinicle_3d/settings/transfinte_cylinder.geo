//+
SetFactory("OpenCASCADE");

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Height = 10;
Radius = 6;

Nb = 20+1;
Rb = 1;
Nc1 = 10+1;
Rc1 = 1.00;
Nz = 2;

//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Radius*0.70710678, Radius*0.70710678, 0, 1.0};
//+
Point(3) = {-Radius*0.70710678, Radius*0.70710678, 0, 1.0};
//+
Point(4) = {-Radius*0.70710678, -Radius*0.70710678, 0, 1.0};
//+
Point(5) = {Radius*0.70710678, -Radius*0.70710678, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Line(5) = {1, 2};
//+
Line(6) = {1, 3};
//+
Line(7) = {1, 4};
//+
Line(8) = {1, 5};
//+
Transfinite Curve {1} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {2} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {3} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {4} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {5} = Nb Using Progression Rb;
//+
Transfinite Curve {6} = Nb Using Progression Rb;
//+
Transfinite Curve {7} = Nb Using Progression Rb;
//+
Transfinite Curve {8} = Nb Using Progression Rb;
//+
Curve Loop(1) = {5, 1, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 2, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 3, -8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 4, -5};
//+
Plane Surface(4) = {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Extrude {0, 0, Height} {
  Surface{4}; Surface{1}; Surface{2}; Surface{3}; 
  Layers{Nz};
  Recombine;
}
//+
Physical Surface("Bottom", 1) = {3, 2, 4, 1};
//+
Physical Surface("Top", 2) = {11, 14, 8, 16};
//+
Physical Surface("Side", 3) = {6, 15, 12, 9};
//+
Physical Volume("Sea", 4) = {3, 4, 1, 2};
