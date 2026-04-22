// Gmsh project created on Tue Mar 17 14:46:21 2026
//SetFactory("OpenCASCADE");
cl__1 = 35;
Point(1) = {30, 8, 0, cl__1};
Point(5) = {25, 15, 0, cl__1};
Point(12) = {0, 0, 0, cl__1};
Point(13) = {0, 77, 0, cl__1};
Point(14) = {66, 77, 0, cl__1};
Point(15) = {66, 0, 0, cl__1};
//
Point(16) = {4, 45, 0, cl__1};
Point(17) = {35, 40, 0, cl__1};
//
Point(18) = {5, 10, 0, cl__1};
Point(19) = {15, 16, 0, cl__1};
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

//
// Fractura 1
Point(101) = {10,25,0,cl__1}; Point(102) = {15,18,0,cl__1};
Line(101) = {101,102}; Line(102) = {102,101};
Curve Loop(1)={101,102};
Plane Surface(1)={1};
// Fractura 2
Point(103) = {20,15,0,cl__1}; Point(104) = {25,22,0,cl__1};
Line(103) = {103,104}; Line(104) = {104,103};
Curve Loop(2)={103,104};
Plane Surface(2)={2};
// Fractura 3
Point(105) = {45,5,0,cl__1}; Point(106) = {35,18,0,cl__1};
Line(105) = {105,106}; Line(106) = {106,105};
Curve Loop(3)={103,104};
Plane Surface(3)={3};
// Fractura 4
Point(107) = {40,27,0,cl__1}; Point(108) = {45,20,0,cl__1};
Line(107) = {107,108}; Line(108) = {108,107};
Curve Loop(4)={107,108};
Plane Surface(4)={4};
// Fractura 5
Point(109) = {50,10,0,cl__1}; Point(110) = {55,18,0,cl__1};
Line(109) = {109,110}; Line(110) = {110,109};
Curve Loop(5)={109,110};
Plane Surface(5)={5};
// Fractura 6
Point(111) = {10,30,0,cl__1}; Point(112) = {15,38,0,cl__1};
Line(111) = {111,112}; Line(112) = {112,111};
Curve Loop(6)={111,112};
Plane Surface(6)={6};

//
//Line(18) = {16,17};
//Line(19) = {17,16};
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
//Curve Loop(12) = {18,19};
//Plane Surface(12) = {12};
//
Curve Loop(13) = {20,21};
Plane Surface(13) = {13};
//
Curve Loop(14) = {22,23};
Plane Surface(14) = {14};

//Curve Loop(1000) = {1, 2, 3, 4, 5, 6, 7, 8, 9, -15, -14, -13, -12};
Curve Loop(1000) = { 16,17,20,21,22,23,101,102,103,104,105,106,107,108,109,110,
111,112,-15, -14, -13, -12};
//Curve Loop(1001) = { 20,21,22,23};

//Curve Loop(1000) = {  -15, -14, -13, -12};

Plane Surface(1000) = {1000};
Physical Surface("k11", 1) = {1000};
//+
Physical Curve("inlet", 2) = {14};
//+
Physical Curve("outlet", 3) = {12};
Physical Curve("noflux", 4) = {13,15};
//Physical Surface("Vugs",6)={11};
Physical Curve("SmallFract",6)={16,17,20,21,22,23,101,102,103,104,105,106,107,108,109,110,
111,112};

Coherence;