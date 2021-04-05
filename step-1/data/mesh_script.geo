// Gmsh project created on Sun Apr 04 11:41:39 2021
Inner_Radius = 5;
Radius = 20;
Height = 20;

Nb = 2; Rb = 1.00;
Nc1 = 33; Rc1 = 1.00;
Nc2 = 33; Rc2 = 1.00;
Nz = 2;

SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {Inner_Radius*0.70710678, Inner_Radius*0.70710678, 0, 1.0};
//+
Point(3) = {-Inner_Radius*0.70710678, Inner_Radius*0.70710678, 0, 1.0};
//+
Point(4) = {Inner_Radius*0.70710678, -Inner_Radius*0.70710678, 0, 1.0};
//+
Point(5) = {-Inner_Radius*0.70710678, -Inner_Radius*0.70710678, 0, 1.0};
//+
Circle(1) = {5, 1, 4};
//+
Circle(2) = {4, 1, 2};
//+
Circle(3) = {2, 1, 3};
//+
Circle(4) = {3, 1, 5};
//+
Point(6) = {Radius*0.70710678, Radius*0.70710678, 0, 1.0};
//+
Point(7) = {-Radius*0.70710678, Radius*0.70710678, 0, 1.0};
//+
Point(8) = {-Radius*0.70710678, -Radius*0.70710678, 0, 1.0};
//+
Point(9) = {Radius*0.70710678, -Radius*0.70710678, 0, 1.0};
//+
Circle(5) = {8, 1, 9};
//+
Circle(6) = {9, 1, 6};
//+
Circle(7) = {6, 1, 7};
//+
Circle(8) = {7, 1, 8};
//+
Line(9) = {5, 8};
//+
Line(10) = {4, 9};
//+
Line(11) = {2, 6};
//+
Line(12) = {3, 7};
//+
Transfinite Curve {1} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {2} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {3} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {4} = Nc1 Using Progression Rc1;
//+
Transfinite Curve {12} = Nb Using Progression Rb;
//+
Transfinite Curve {11} = Nb Using Progression Rb;
//+
Transfinite Curve {10} = Nb Using Progression Rb;
//+
Transfinite Curve {9} = Nb Using Progression Rb;
//+
Transfinite Curve {8} = Nc2 Using Progression Rc2;
//+
Transfinite Curve {7} = Nc2 Using Progression Rc2;
//+
Transfinite Curve {6} = Nc2 Using Progression Rc2;
//+
Transfinite Curve {5} = Nc2 Using Progression Rc2;
//+
Curve Loop(1) = {9, 5, -10, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 6, -11, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 7, -12, -3};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 8, -9, -4};
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
Physical Surface("Top", 1) = {20, 9, 17, 13};
//+
Physical Surface("Bottom", 2) = {3, 2, 4, 1};
//+
Physical Surface("Side", 3) = {18, 14, 10, 6};
//+
Physical Surface("Inner_Side", 4) = {16, 19, 8, 12};
//+
Physical Volume("Sea", 5) = {4, 3, 2, 1};
