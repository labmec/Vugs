cl__1 = 35;
Point(1) = {1.25, .25, 0, cl__1};
Point(5) = {1.25, 0.75, 0, cl__1};
Point(12) = {0.5, 0, 0, cl__1};
Point(13) = {0.5, 1, 0, cl__1};
Point(14) = {2, 1, 0, cl__1};
Point(15) = {2, 0, 0, cl__1};


Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 12};
Line(16) = {5,1};
Line(17) = {1,5};
//Curve Loop(11) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//Plane Surface(11) = {11};
Curve Loop(11) = {16,17};
Plane Surface(11) = {11};
Physical Curve("SmallFract",5)={16,17};

//Curve Loop(1000) = {1, 2, 3, 4, 5, 6, 7, 8, 9, -15, -14, -13, -12};
Curve Loop(1000) = { 16,17, -15, -14, -13, -12};
//Curve Loop(1000) = {  -15, -14, -13, -12};

Plane Surface(1000) = {1000};
Physical Surface("k11", 1) = {1000};
//+
Physical Curve("inlet", 2) = {14};
//Physical Point("FracEnds",10)={1,5};

