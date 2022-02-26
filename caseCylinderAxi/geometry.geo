// Gmsh project created on Thu Feb 24 08:19:53 2022
//+
Point(1) = {0, 0, 0, 0.025};
//+
Point(2) = {-0.5, 0, 0, 0.025};
//+
Point(3) = {0, 0.5, 0, 0.025};
//+
Point(4) = {0, 2, 0, 0.025};
//+
Point(5) = {-1, 1, 0, 0.025};
//+
Point(6) = {-1, 0, 0, 0.025};
//+
Circle(1) = {2, 1, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 2};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("wall", 6) = {1};
//+
Physical Curve("outlet", 7) = {2};
//+
Physical Curve("inlet", 8) = {3, 4};
//+
Physical Curve("sym", 9) = {5};
