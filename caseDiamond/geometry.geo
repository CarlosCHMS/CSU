// Gmsh project created on Tue Feb 22 15:29:44 2022
//+
Point(1) = {-2, -2, 0, 0.1};
//+
Point(2) = {2, -2, 0, 0.1};
//+
Point(3) = {2, 2, 0, 0.1};
//+
Point(4) = {-2, 2, 0, 0.1};
//+
Point(5) = {-0.5, 0, 0, 0.05};
//+
Point(6) = {0, 0.05, 0, 0.05};
//+
Point(7) = {0.5, 0, 0, 0.05};
//+
Point(8) = {0, -0.05, 0, 0.05};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("inlet") = {4};
//+
Physical Curve("outlet") = {2};
//+
Physical Curve("wallFar") = {3, 1};
//+
Physical Curve("wall") = {8, 5, 6, 7};
