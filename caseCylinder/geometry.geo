// Gmsh project created on Thu Feb 24 08:19:53 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.05};
//+
Point(2) = {-1, 0, 0, 0.05};
//+
Point(3) = {0, 1, 0, 0.05};
//+
Point(4) = {0, 3, 0, 0.5};
//+
Point(5) = {-3, 0, 0, 0.5};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {4, 1, 5};
//+
Line(3) = {5, 2};
//+
Line(4) = {3, 4};
//+
Physical Curve("wall", 5) = {1};
//+
Physical Curve("outlet", 6) = {4};
//+
Physical Curve("inlet", 7) = {2};
//+
Physical Curve("sym", 8) = {3};
//+
Curve Loop(1) = {3, 1, 4, 2};
//+
Plane Surface(1) = {1};
