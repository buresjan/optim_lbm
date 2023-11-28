Merge "DEFINE_STL_FILE";
//+
SetFactory("OpenCASCADE");

X = DEFINE_X;
Y = DEFINE_Y;
Z = DEFINE_Z;
DX = DEFINE_DX;
DY = DEFINE_DY;
DZ = DEFINE_DZ;
VOXEL_SIZE = DEFINE_VOXEL_SIZE;

Box(1) = {X - VOXEL_SIZE, Y - VOXEL_SIZE, Z - VOXEL_SIZE, DX + 2*VOXEL_SIZE, DY + 2*VOXEL_SIZE, DZ + 2*VOXEL_SIZE};