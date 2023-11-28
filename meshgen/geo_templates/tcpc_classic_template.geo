SetFactory("OpenCASCADE");

h = DEFINE_H;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

//RES = DEFINE_RESOLUTION;
RES = 1.0;
ANGLE = DEFINE_ANGLE;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Model of extracardiac conduit /////////////////////////////////////////////
/////////////////////////////////////////////// ID - 00001 //////////////////////////////////////////////////////

// Points that determine the central axis of the conduit - used for forming a closed loop
Point(1) = {RES * (-1.00), RES * (0.64 + 0.0000), 0.00};
Point(9) = {RES * ( 0.25), RES * (0.64 + 0.0000), 0.00};

// Sub-level point for correction during rotation of the conduit
Point(2) = {RES * (-1.00), RES * (0.64 + 0.0108), 0.0};

// Main points forming the spline to be rotated and revolved
Point(3) = {RES * ( 0.00), RES * (0.64 + 0.0408), 0.0};
Point(4) = {RES * ( 0.05), RES * (0.64 + 0.0440), 0.0};
Point(5) = {RES * ( 0.10), RES * (0.64 + 0.0475), 0.0};
Point(6) = {RES * ( 0.15), RES * (0.64 + 0.0490), 0.0};
Point(7) = {RES * ( 0.20), RES * (0.64 + 0.0530), 0.0};
Point(8) = {RES * ( 0.25), RES * (0.64 + 0.0780), 0.0};

// Lines that form a close loop
Line(1) = {8, 9};
Line(2) = {9, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

// Spline connecting the main points
Spline(5) = {3, 4, 5, 6, 7, 8};

// Closed loop formed by the lines
Curve Loop(1) = {1, 2, 3, 4, 5};

// Make a surface inside of the closed loop
Plane Surface(1) = {1};

// Revolve the surface around the axis X that is translated to the central point of the conduit creating Volume{1}
Extrude { {1, 0, 0}, {RES * (-1.00), RES * (0.64 + 0.0000), 0.00} , 2*Pi } {
  Surface{1}; Recombine;
}

// Rotate the created volume by specified angle around the axis Z that is translated to the central point of the conduit
Rotate { {0, 0, 1}, {RES * (-1.00), RES * (0.64 + 0.0000), 0.00}, ANGLE} {
  Volume{1};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Model of pulmonary artery /////////////////////////////////////////////
/////////////////////////////////////////////// ID - 00100 /////////////////////////////////////////////////////

// Points that determine the central axis of the pulmonary artery - used for forming a closed loop
Point(101) = {RES * (0.30),         0.00, 0.00};
Point(108) = {RES * (0.30), RES * (1.28), 0.00};

// Main points forming the spline to be revolved
Point(102) = {RES * (0.30 + 0.090),         0.000, 0.00};
Point(103) = {RES * (0.30 + 0.095), RES * (0.256), 0.00};
Point(104) = {RES * (0.30 + 0.100), RES * (0.512), 0.00};
Point(105) = {RES * (0.30 + 0.090), RES * (0.768), 0.00};
Point(106) = {RES * (0.30 + 0.100), RES * (1.024), 0.00};
Point(107) = {RES * (0.30 + 0.085), RES * (1.280), 0.00};

// Lines that form a close loop
Line(101) = {107, 108};
Line(102) = {108, 101};
Line(103) = {101, 102};

// Spline connecting the main points
Spline(104) = {102, 103, 104, 105, 106, 107};

// Closed loop formed by the lines
Curve Loop(101) = {101, 102, 103, 104};

// Make a surface inside of the closed loop
Plane Surface(101) = {101};

// Revolve the surface around the axis Y that is translated to the central point of the pulmonary artery creating Volume{2}
Extrude { {0, 1, 0} , {RES * (0.30), 0.0, 0.0} , 2*Pi } {
  Surface{101}; Recombine;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// Model of vena cava superior ///////////////////////////////////////////
/////////////////////////////////////////////// ID - 10000 /////////////////////////////////////////////////////

// Points that determine the central axis of the vena cava - used for forming a closed loop
Point(10001) = {RES * (0.64), RES * (0.64), 0.0};
Point(10007) = {RES * (0.36), RES * (0.64), 0.0};

// Main points forming the spline to be revolved
Point(10002) = {RES * (0.64), RES * (0.64 + 0.0408), 0.0};
Point(10003) = {RES * (0.57), RES * (0.64 + 0.0440), 0.0};
Point(10004) = {RES * (0.50), RES * (0.64 + 0.0475), 0.0};
Point(10005) = {RES * (0.43), RES * (0.64 + 0.0490), 0.0};
Point(10006) = {RES * (0.36), RES * (0.64 + 0.0530), 0.0};

// Lines that form a close loop
Line(10001) = {10002, 10001};
Line(10002) = {10001, 10007};
Line(10003) = {10007, 10006};

// Spline connecting the main points
Spline(10004) = {10002, 10003, 10004, 10005, 10006};

// Closed loop formed by the lines
Curve Loop(10001) = {10001, 10002, 10003, 10004};

// Make a surface inside of the closed loop
Plane Surface(10001) = {10001};

// Revolve the surface around the axis X that is translated to the central point of the vena cava creating Volume{3}
Extrude { {1, 0, 0} , {0, RES * (0.64), 0} , 2*Pi } {
  Surface{10001}; Recombine;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Cleaning up and merging /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Merge the three volumes creating one Volume {4} while deleting the three partial volumes
BooleanUnion(4) = { Volume{1}; Delete; }{ Volume{2, 3}; Delete; };

// Delete everything that is left -->
Recursive Delete {
  Surface{1, 101, 10001};
}

Recursive Delete {
  Point{4, 5, 6, 7, 102, 103, 104, 105, 106, 107, 10002, 10003, 10004, 10005};
}
// Delete everything that is left <--

// Create a box domain that is used to slice potential "tilted" bottom part of the rotated conduit - we want it to be
//  parallel with the ZY plane
Box(5) = {0.0, 0.0, RES * (-0.101), RES * 0.64, RES * 1.28, 2 * RES * 0.101};

// Intersect the box with the Volume {4}
BooleanIntersection{ Volume{5}; Delete;}{ Volume{4}; Delete;};