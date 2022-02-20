// Gmsh project created on Sun Feb 20 15:04:05 2022
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wall1") = {1};
//+
Physical Curve("wall2") = {2};
//+
Physical Curve("wall3") = {3};
//+
Physical Curve("wall4") = {4};
//+
Physical Surface("surf") = {1};
