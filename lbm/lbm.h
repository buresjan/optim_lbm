#pragma once

#include "defs.h"
#include "lbm_data.h"

template< 
	typename LBM_TYPE,
	typename MACRO,
	typename CPU_MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
struct LBM
{
	using T_TRAITS = typename LBM_TYPE::T_TRAITS;
	using idx = typename T_TRAITS::idx;
	using dreal = typename T_TRAITS::dreal;
	using real = typename T_TRAITS::real;
	using map_t = typename T_TRAITS::map_t;

	map_t null=0;
	dreal drealnull=-1e88;

	// KernelData contains only the necessary data for the CUDA kernel. these are copied just before the kernel is called
	LBM_DATA data;
	// block size for the LBM kernel (used in core.h)
	int block_size;

	map_t &map(idx gi) { if (gi>=0 && gi < X*Y*Z) return hmap[gi]; else return null; }
	map_t &map(idx x, idx y, idx z) { return map(pos(x,y,z)); }

	dreal *cpumacro;
	dreal *dmacro, *hmacro;
	map_t *hmap, *dmap;
	dreal *dfs[DFMAX], *hfs[DFMAX];
//	dreal *lat; // pointer to hfs[0] for backward compatibility

	// macroscopic quantities in RAM
	bool *wall; // indicates whether there is a wall (custom define outside the class State)

	idx pos(idx x, idx y, idx z) { return POS(x, y, z, X, Y, Z); }
	idx mpos(int variable_index, idx x, idx y, idx z) { return mPOS(variable_index, x, y, z, X, Y, Z); }

	dreal& macro(int index, idx x, idx y, idx z) { if (index>=MACRO::N) return drealnull; else return hmacro[mpos(index,x,y,z)]; }
	dreal& cpu_macro(int index, idx x, idx y, idx z) { if (index>=CPU_MACRO::N) return drealnull; else return cpumacro[mpos(index,x,y,z)]; }


	// LBM parameters (*copied from StateData)
	idx X, Y, Z;
//	size_t size_dreal, size_real;

	// parameters for multi-GPU
	bool use_multiple_gpus = true;
	idx offset_X = 0;
	idx local_X = 0;
	idx overlap_left = 0;
	idx overlap_right = 0;

	// input parameters: constant in time
	real physDl; 				// spatial step (fixed throughout the computation)
	real physDt;				// temporal step (fixed or variable throughout the computation)
	real physFilDl;				// spatial step for filaments(should be fixed throught the computation)
	real physViscosity;			// physical viscosity of the fluid
	real physFluidDensity;			// physical (characteristic) density of the fluid (constant)
	real physFinalTime;			// default 1e10
	real physCharLength;			// characteristic length, default is physDl * (real)Y but you can specify that manually
	int iterations;			// number of lbm iterations

	bool terminate;			// flag for terminal error detection

//	bool noFlow;			// trigers no flow boundary condition

//	real Re(real physvel) { return fabs(physvel) * physDl * (real)Y / physViscosity; } // TODO: change Y to charLength --- specify this explicitely
	real Re(real physvel) { return fabs(physvel) * physCharLength / physViscosity; } // TODO: change Y to charLength --- specify this explicitely
	dreal lbmViscosity() { return (dreal) (physDt / physDl / physDl * physViscosity); }
	real physTime() { return physDt*(real)iterations; }
	real lbm2physVelocity(real lbm_velocity) { return lbm_velocity / physDt * physDl; }
	real lbm2physForce(real lbm_force) { return lbm_force * physDl / physDt / physDt; }
	real phys2lbmVelocity(real phys_velocity)  { return phys_velocity * physDt / physDl; }
	real phys2lbmForce(real phys_force) { return phys_force / physDl * physDt * physDt; }
	dreal lbmInputDensity;
	
//	real physNormVelocity(idx gi) { return NORM( lbm2physVelocity(hvx[gi]), lbm2physVelocity(hvy[gi]), lbm2physVelocity(hvz[gi]) ); }
//	real physDensity(idx gi) { return hrho[gi]*physFluidDensity; }
//	real physNormVelocity(idx GX, idx GY, idx GZ) { return physNormVelocity(pos(GX,GY,GZ)); }
//	real physDensity(idx GX, idx GY, idx GZ) { return physDensity(pos(GX,GY,GZ)); }

	void resetForces() { resetForces(0,0,0);}
	void resetForces(real ifx, real ify, real ifz);
	void copyForcesToDevice();
	
	// all this is needed for IBM only and forcing
	dreal* hrho() { return MACRO::rho(hmacro, X*Y*Z); }
	dreal* hvx() { return MACRO::vx(hmacro, X*Y*Z); }
	dreal* hvy() { return MACRO::vy(hmacro, X*Y*Z); }
	dreal* hvz() { return MACRO::vz(hmacro, X*Y*Z); }
	dreal* hfx() { return MACRO::fx(hmacro, X*Y*Z); }
	dreal* hfy() { return MACRO::fy(hmacro, X*Y*Z); }
	dreal* hfz() { return MACRO::fz(hmacro, X*Y*Z); }

	dreal& hrho(idx i) { return hrho()[i]; }
	dreal& hvx(idx i)  { return hvx()[i]; }
	dreal& hvy(idx i)  { return hvy()[i]; }
	dreal& hvz(idx i)  { return hvz()[i]; }
	dreal& hfx(idx i)  { return hfx()[i]; }
	dreal& hfy(idx i)  { return hfy()[i]; }
	dreal& hfz(idx i)  { return hfz()[i]; }

	dreal* drho() { return MACRO::rho(dmacro, X*Y*Z); }
	dreal* dvx() { return MACRO::vx(dmacro, X*Y*Z); }
	dreal* dvy() { return MACRO::vy(dmacro, X*Y*Z); }
	dreal* dvz() { return MACRO::vz(dmacro, X*Y*Z); }
	dreal* dfx() { return MACRO::fx(dmacro, X*Y*Z); }
	dreal* dfy() { return MACRO::fy(dmacro, X*Y*Z); }
	dreal* dfz() { return MACRO::fz(dmacro, X*Y*Z); }

	void copyMapToDevice();
	void copyMapToHost();
	void copyMacroToHost();
	void copyMacroToDevice();
	void copyDFsToHost(dreal*source, dreal*target);
	void copyDFsToDevice(dreal*source, dreal*target);
	void copyDFsToHost();
	void copyDFsToDevice();
//	void copyEqLatticeToDevice();

	void computeCPUMacroFromLat();

	void clearHostMacro();
	void synchronizeDFsToRight(int gpu_id, int gpu_id_right, LBM& lbm_right);
	void synchronizeDFsToLeft(int gpu_id, int gpu_id_left, LBM& lbm_left);
	void synchronizeMacroToLeft(int gpu_id, int gpu_id_left, LBM& lbm_left);
	void synchronizeMacroToRight(int gpu_id, int gpu_id_right, LBM& lbm_right);

	void defineWall(idx x, idx y, idx z, bool value);
	void projectWall() { for (idx i=0; i < X*Y*Z; i++) if (wall[i]) map(i) = LBM_BC::GEO_WALL; };
	void setBoundaryX(idx x, map_t value) { for (idx y=0;y<Y;y++) for (idx z=0;z<Z;z++) map(x,y,z)=value; }
	void setBoundaryY(idx y, map_t value) { for (idx x=0;x<X;x++) for (idx z=0;z<Z;z++) map(x,y,z)=value; }
	void setBoundaryZ(idx z, map_t value) { for (idx x=0;x<X;x++) for (idx y=0;y<Y;y++) map(x,y,z)=value; }

	inline bool getWall(idx x, idx y, idx z)
	{
		if (x<0 || x>=X || y <0 || y>=Y || z<0 || z>=Z) return false;
		return wall[pos(x,y,z)];
	}
	inline bool isFluid(idx x, idx y, idx z)
	{
		if (x<0 || x >= X || y<0 || y >= Y || z<0 || z>=Z) return false;
		return LBM_BC::isFluid(map(x,y,z));
	}
	void resetMap(map_t geo_type);
	void resetVelocities(real ivx, real ivy, real ivz);
	void resetLattice(real irho, real ivx, real ivy, real ivz);

	bool quit() { return terminate; }

	void allocateDeviceData();
	void updateKernelData(); 		// copy physical parameters to data structure accessible by the CUDA kernel

	LBM(idx iX, idx iY, idx iZ, real iphysViscosity, real iphysDl, real iphysDt);
	~LBM();
};

#include "lbm.hpp"
