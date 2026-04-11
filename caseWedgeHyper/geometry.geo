// Gmsh project created on Tue Feb 22 15:29:44 2022
l = 0.0075;
//+
Point(1) = {-0.2, 0, 0, l};
//+
Point(2) = {0, 0, 0, l};
//+
Point(3) = {0.8, 0.291176, 0, l};
//+
Point(4) = {0.8, 0.7, 0, l};
//+
Point(5) = {-0.2, 0.3, 0, l};
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
Physical Curve("wall") = {1, 2};
//+
Physical Curve("out") = {3};
//+
Physical Curve("in1") = {4};
//+
Physical Curve("in2") = {5};
