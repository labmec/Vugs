//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {70, 0, 0, 1.0};
//+
Point(3) = {70, 80, 0, 1.0};
//+
Point(4) = {0, 80, 0, 1.0};
//+
Point(5) = {20, 40, 0, 1.0};
//+
Point(6) = {50, 40, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Point("FracEnds", 10) = {5, 6};
//+
Physical Curve("outlet", 3) = {2};
//+
Physical Curve("inlet", 2) = {4};
//+
Physical Curve("noflux", 4) = {3, 1};
//+
Physical Surface("k11", 1) = {1};
//+
Physical Curve("SmallFract", 5) = {5,6};

// Define number of divisions on opposite curves
Transfinite Curve {1, 3} = 10; 
Transfinite Curve {2, 4} = 15;

// Apply transfinite constraint to the surface
Transfinite Surface {1};

// Convert the default triangles into quads
Recombine Surface {1};
