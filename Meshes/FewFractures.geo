// Gmsh project created on Tue Mar 17 14:46:21 2026
//SetFactory("OpenCASCADE");
cl__1 = 10;
Point(1) = {30, 8, 0, cl__1};
Point(5) = {25, 15, 0, cl__1};
Point(12) = {0, 0, 0, cl__1};
Point(13) = {0, 77, 0, cl__1};
Point(14) = {66, 77, 0, cl__1};
Point(15) = {66, 0, 0, cl__1};
//
Point(16) = {40, 25, 0, cl__1};
Point(17) = {35, 40, 0, cl__1};
//
Point(18) = {5, 10, 0, cl__1};
Point(19) = {15, 21, 0, cl__1};
//
Point(20) = {55, 60, 0, cl__1};
Point(21) = {45, 68, 0, cl__1};
//
//Point(22) = {66, 0, 0, cl__1};
//Point(23) = {66, 0, 0, cl__1};

Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 12};
Line(16) = {5,1};
Line(17) = {1,5};


Line(18) = {16,17};
Line(19) = {17,16};
//
Line(20) = {18,19};
Line(21) = {19,18};
//
Line(22) = {20,21};
Line(23) = {21,20};

//Curve Loop(11) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//Plane Surface(11) = {11};
//
//
Curve Loop(11) = {16,17};
Plane Surface(11) = {11};
//
//
//
Curve Loop(12) = {18,19};
Plane Surface(12) = {12};
//
Curve Loop(13) = {20,21};
Plane Surface(13) = {13};
//
Curve Loop(14) = {22,23};
Plane Surface(14) = {14};

//Curve Loop(1000) = {1, 2, 3, 4, 5, 6, 7, 8, 9, -15, -14, -13, -12};
Curve Loop(1000) = { 16,17,18,19,20,21,22,23,-15, -14, -13, -12};
//Curve Loop(1001) = { 20,21,22,23};

//Curve Loop(1000) = {  -15, -14, -13, -12};

Plane Surface(1000) = {1000};
Physical Surface("k11", 1) = {1000};
//+
Physical Curve("inlet", 2) = {14};
//+
Physical Curve("outlet", 3) = {12};
Physical Curve("top", 4) = {13,15};
//Physical Curve("bottom", 7) = {15};
//Physical Surface("Vugs",6)={11};
Physical Curve("SmallFract",5)={16,17,18,19,20,21,22,23};
Physical Point("FracEnds",10)={1,5,16,17,18,19,20,21};

Coherence;

//+
Physical Curve(" bottom", 5) -= {15, 13};
