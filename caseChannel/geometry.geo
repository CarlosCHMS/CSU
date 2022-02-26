// Gmsh project created on Sat Feb 26 11:46:45 2022
//+
Point(1) = {0, -2.958, 0, 0.05};
//+
Point(2) = {0.5, 0, 0, 0.05};
//+
Point(3) = {1.5, 0, 0, 0.05};
//+
Point(4) = {1.5, 1, 0, 0.05};
//+
Point(5) = {-1.5, 1, 0, 0.05};
//+
Point(6) = {-1.5, 0, 0, 0.05};
//+
Point(7) = {-0.5, 0, 0, 0.05};
//+
Circle(1) = {7, 1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet") = {5};
//+
Physical Curve("outlet") = {3};
//+
Physical Curve("wallDown") = {6, 1, 2};
//+
Physical Curve("wallUp") = {4};
