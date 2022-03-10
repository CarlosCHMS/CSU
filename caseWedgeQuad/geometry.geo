// Gmsh project created on Tue Feb 22 15:29:44 2022
//+
Point(1) = {-0.2, 0, 0, 0.05};
//+
Point(2) = {0, 0, 0, 0.05};
//+
Point(3) = {0.8, 0.291176, 0, 0.05};
//+
Point(4) = {0.8, 1.5, 0, 0.05};
//+
Point(5) = {-0.2, 1.5, 0, 0.05};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Curve Loop(1) = {2, 3, 4, 5, 1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wall1") = {1, 2};
//+
Physical Curve("wall2") = {3};
//+
Physical Curve("wall3") = {4};
//+
Physical Curve("wall4") = {5};
