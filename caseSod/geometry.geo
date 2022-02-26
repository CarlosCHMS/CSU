// Gmsh project created on Sat Feb 26 16:54:21 2022
//+
Point(1) = {0, 0, 0, 0.4};
//+
Point(2) = {50, 0, 0, 0.4};
//+
Point(3) = {50, 5, 0, 0.4};
//+
Point(4) = {0, 5, 0, 0.4};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wallDown") = {1};
//+
Physical Curve("wallRight") = {2};
//+
Physical Curve("wallUp") = {3};
//+
Physical Curve("wallLeft") = {4};
