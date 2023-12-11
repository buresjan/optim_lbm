SetFactory("OpenCASCADE");

h = DEFINE_H;
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

OFFSET = DEFINE_OFFSET;

// Middle rectangle
Rectangle(1) = {0, 0.20, 0, 1.00, 0.10, 0};

// Lower rectangle
Rectangle(2) = {0.45 + OFFSET, 0, 0, 0.10, 0.20, 0};

// Upper rectangle
Rectangle(3) = {0.45, 0.30, 0, 0.10, 0.20, 0};

BooleanUnion(4) = { Surface{1}; Delete; }{ Surface{2, 3}; Delete; };
