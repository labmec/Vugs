// =====================================
// PARÁMETROS
// =====================================
cl__1 = 300;

// =====================================
// NODOS DEL CUBO
// =====================================
Point(1) = {0, 0, 0, cl__1};
Point(2) = {616, 0, 0, cl__1};
Point(3) = {616, 676, 0, cl__1};
Point(4) = {0, 676, 0, cl__1};

Point(5) = {0, 0, 910, cl__1};
Point(6) = {616, 0, 910, cl__1};
Point(7) = {616, 676, 910, cl__1};
Point(8) = {0, 676, 910, cl__1};

// =====================================
// ARISTAS
// =====================================
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// =====================================
// CARAS DEL CUBO
// =====================================
Line Loop(13) = {1,2,3,4};
Plane Surface(14) = {13};

Line Loop(15) = {5,6,7,8};
Plane Surface(16) = {15};

Line Loop(17) = {1,10,-5,-9};
Plane Surface(18) = {17};

Line Loop(19) = {2,11,-6,-10};
Plane Surface(20) = {19};

Line Loop(21) = {3,12,-7,-11};
Plane Surface(22) = {21};

Line Loop(23) = {4,9,-8,-12};
Plane Surface(24) = {23};

// =====================================
// VOLUMEN
// =====================================
Surface Loop(30) = {14,16,18,20,22,24};
Volume(31) = {30};

// =====================================
// FRACTURA (PLANO INTERNO)
// =====================================
// Definimos 4 puntos dentro del cubo
Point(100) = {200, 200, 450, cl__1};
Point(101) = {400, 200, 350, cl__1};
Point(102) = {400, 400, 350, cl__1};
Point(103) = {200, 400, 450, cl__1};

// Líneas de la fractura
Line(100) = {100,101};
Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,100};

// Superficie de la fractura
Line Loop(104) = {100,101,102,103};
Plane Surface(105) = {104};

// =====================================
// EMBEBER FRACTURA EN EL VOLUMEN
// =====================================
Surface{105} In Volume{31};

// =====================================
// PHYSICAL GROUPS
// =====================================
Physical Volume("CuboExterno", 1) = {31};

Physical Surface("Fracture", 2) = {105};

// Boundaries (ajusta según tu solver)
Physical Surface("inlet", 3) = {20};
Physical Surface("outlet", 4) = {24};
Physical Surface("noflux", 5) = {14,16,18,22};

// =====================================
//Mesh 3;