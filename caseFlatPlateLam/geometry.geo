// Gmsh project created on Sat Feb 26 16:54:21 2022
//+
Point(1) = {0, 0, 0, 0.4};
//+
Point(2) = {0.5, 0, 0, 0.4};
//+
Point(3) = {0.5, 0.006, 0, 0.4};
//+
Point(4) = {0, 0.006, 0, 0.4};
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
Physical Curve("wall") = {1};
//+
Physical Curve("outlet") = {2, 3};
//+
Physical Curve("inlet") = {4};
//+
Transfinite Curve{1} = 200 Using Progression 1.02;
Transfinite Curve{3} = 200 Using Progression 0.98039;

//Transfinite Curve{2} = 40; 
//Transfinite Curve{4} = 40;

Transfinite Curve{2} = 15 Using Progression 1.1;
Transfinite Curve{4} = 15 Using Progression 0.909090;

Transfinite Surface{1};
