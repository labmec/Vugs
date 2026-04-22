// Gmsh project
cl__1 = 35;

// Dominio
Point(12) = {0, 0, 0, cl__1};
Point(13) = {0, 77, 0, cl__1};
Point(14) = {66, 77, 0, cl__1};
Point(15) = {66, 0, 0, cl__1};

// Contorno
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 12};

// ==========================
// 15 FRACTURAS (pares de líneas)
// ==========================

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
Point(105) = {30,5,0,cl__1}; Point(106) = {35,18,0,cl__1};
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
// Fractura 7
Point(113) = {10,45,0,cl__1}; Point(114) = {25,42,0,cl__1};
Line(113) = {113,114}; Line(114) = {114,113};
Curve Loop(7)={113,114};
Plane Surface(7)={7};
// Fractura 8
Point(115) = {30,32,0,cl__1}; Point(116) = {35,40,0,cl__1};
Line(115) = {115,116}; Line(116) = {116,115};
Curve Loop(8)={115,116};
Plane Surface(8)={8};
// Fractura 9
Point(117) = {40,40,0,cl__1}; Point(118) = {50,38,0,cl__1};
Line(117) = {117,118}; Line(118) = {118,117};
Curve Loop(9)={117,118};
Plane Surface(9)={9};
// Fractura 10
Point(119) = {50,35,0,cl__1}; Point(120) = {55,42,0,cl__1};
Line(119) = {119,120}; Line(120) = {120,119};
Curve Loop(10)={119,120};
Plane Surface(10)={10};
// Fractura 11
Point(121) = {3,65,0,cl__1}; Point(122) = {15,65,0,cl__1};
Line(121) = {121,122}; Line(122) = {122,121};
Curve Loop(11)={121,122};
Plane Surface(11)={11};
// Fractura 12
Point(123) = {14,60,0,cl__1}; Point(124) = {35,70,0,cl__1};
Line(123) = {123,124}; Line(124) = {124,123};
Curve Loop(12)={123,124};
Plane Surface(12)={12};
// Fractura 13
Point(125) = {30,50,0,cl__1}; Point(126) = {40,68,0,cl__1};
Line(125) = {125,126}; Line(126) = {126,125};
Curve Loop(13)={125,126};
Plane Surface(13)={13};
// Fractura 14
Point(127) = {48,44,0,cl__1}; Point(128) = {45,65,0,cl__1};
Line(127) = {127,128}; Line(128) = {128,127};
Curve Loop(14)={127,128};
Plane Surface(14)={14};
// Fractura 15
Point(129) = {50,60,0,cl__1}; Point(130) = {55,70,0,cl__1};
Line(129) = {129,130}; Line(130) = {130,129};
Curve Loop(15)={129,130};
Plane Surface(15)={15};
// ==========================
// SUPERFICIE PRINCIPAL
// ==========================
Curve Loop(1000) = {
101,102,103,104,105,106,107,108,109,110,
111,112,113,114,115,116,117,118,119,120,
121,122,123,124,125,126,127,128,129,130,
-15,-14,-13,-12
};

Plane Surface(1000) = {1000};

// ==========================
// PHYSICAL GROUPS
// ==========================
Physical Surface("k11", 1) = {1000};

Physical Curve("inlet", 2) = {14};
Physical Curve("outlet", 3) = {12};
Physical Curve("noflux", 4) = {13,15};

Physical Curve("SmallFract", 6) = {
101,102,103,104,105,106,107,108,109,110,
111,112,113,114,115,116,117,118,119,120,
121,122,123,124,125,126,127,128,129,130
};

Coherence;