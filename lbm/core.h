#pragma once

#include <omp.h>

#include "defs.h"
#include "state.h"
#include "lbm.h"

// default
#include "lbm_macro.h"



#include "lbm_eq.h"
#include "lbm_eq_inv_cum.h"
#include "lbm_eq_well.h"
#include "lbm_eq_entropic.h"

#include "lbm_streaming.h"

#include "lbm_bc.h"

#include "lbm_col_cum.h"
#include "lbm_col_bgk.h"
#include "lbm_col_clbm.h"
#include "lbm_col_fclbm.h"
#include "lbm_col_mrt.h"
#include "lbm_col_srt.h"
#include "lbm_col_cum_sgs.h"
#include "lbm_col_kbc_n.h"
#include "lbm_col_kbc_c.h"
#include "lbm_col_srt_modif_force.h"
#include "lbm_col_clbm_fei.h"

#include "lbm_col_srt_well.h"
#include "lbm_col_clbm_well.h"
#include "lbm_col_cum_well.h"
#include "lbm_col_bgk_well.h"


template <
	typename LBM_TYPE,
	typename STREAMING,
	typename MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
#ifdef USE_CUDA
__global__ void cudaLBMKernel(LBM_DATA SD)
#else
void LBMKernel(
	LBM_DATA SD, 
	typename LBM_TYPE::T_TRAITS::idx x,
	typename LBM_TYPE::T_TRAITS::idx y,
	typename LBM_TYPE::T_TRAITS::idx z
)
#endif
{
	using dreal = typename LBM_TYPE::T_TRAITS::dreal;
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	using map_t = typename LBM_TYPE::T_TRAITS::map_t;

	#ifdef USE_CUDA
	idx x = threadIdx.x + blockIdx.x * blockDim.x;
	idx y = threadIdx.y + blockIdx.y * blockDim.y;
	idx z = threadIdx.z + blockIdx.z * blockDim.z;
	#endif
	idx gi = POS(x, y, z, SD.X, SD.Y, SD.Z);
	map_t gi_map = SD.dmap[gi];
//	idx gi_macro = POS(x + SD.overlap_left,y,z, SD.X + SD.overlap_left + SD.overlap_right, SD.Y, SD.Z); // FIXME
	idx gi_macro = POS(x + SD.macro_overlap_left,y,z, SD.X + SD.macro_overlap_left + SD.macro_overlap_right, SD.Y, SD.Z);

	KernelStruct<dreal> KS;

	// copy quantities
	MACRO::copyQuantities(SD,KS,gi); // FIXME possible error: + overlap? 

	idx xp,xm,yp,ym,zp,zm;
	if (LBM_BC::isPeriodic(gi_map))
	{
		// handle overlaps between GPUs
		xp = (SD.overlap_right==0 && x == SD.X-1) ? 0 : (x+1);
		xm = (SD.overlap_left==0 && x == 0) ? (SD.X-1) : (x-1);
		yp = (y == SD.Y-1) ? 0 : (y+1);
		ym = (y == 0) ? (SD.Y-1) : (y-1);
		zp = (z == SD.Z-1) ? 0 : (z+1);
		zm = (z == 0) ? (SD.Z-1) : (z-1);
	} else {
		// handle overlaps between GPUs
		xp = (SD.overlap_right>0) ? x+1 : MIN(x+1, SD.X-1);
		xm = (SD.overlap_left>0) ? x-1 : MAX(x-1,0);
		yp = MIN(y+1, SD.Y-1);
		ym = MAX(y-1,0);
		zp = MIN(z+1, SD.Z-1);
		zm = MAX(z-1,0);
	}

	// Streaming
	if (LBM_BC::isStreaming(gi_map))
		STREAMING::streaming(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);

	// compute Density & Velocity
	if (LBM_BC::isComputeDensityAndVelocity(gi_map))
		LBM_TYPE::computeDensityAndVelocity(KS);


	// boundary conditions
	if (LBM_BC::template BC<LBM_TYPE,STREAMING,LBM_DATA>(SD,KS,gi_map,xm,x,xp,ym,y,yp,zm,z,zp)==false)
	{
		LBM_TYPE::collision(KS);
	}

	LBM_TYPE::copyKS2DFout(SD,KS,x,y,z);
//	MACRO::outputMacro(SD.dmacro,SD.X*SD.Y*SD.Z,KS,gi);
	MACRO::outputMacro(SD.dmacro,(SD.X+SD.macro_overlap_left + SD.macro_overlap_right)*SD.Y*SD.Z,KS,gi_macro);
}


// wrapper: work on Macro before LBMKernel
template <
	typename LBM_TYPE,
	typename STREAMING,
	typename MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
#ifdef USE_CUDA
__global__ void cudaMacroWorker(LBM_DATA SD)
#else
void MacroWorker(
	LBM_DATA SD, 
	typename LBM_TYPE::T_TRAITS::idx x,
	typename LBM_TYPE::T_TRAITS::idx y,
	typename LBM_TYPE::T_TRAITS::idx z
)
#endif
{
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	#ifdef USE_CUDA
	idx x = threadIdx.x + blockIdx.x * blockDim.x;
	idx y = threadIdx.y + blockIdx.y * blockDim.y;
	idx z = threadIdx.z + blockIdx.z * blockDim.z;
	#endif
	MACRO::template kernelWorker<LBM_TYPE, STREAMING, LBM_DATA, LBM_BC>(SD,x,y,z);
}

// initial condition --> hmacro on CPU
template <
	typename LBM_TYPE,
	typename LBM_BC,
	typename MACRO,
	typename LBM
>
void LBMKernelInit(
	LBM& lbm,
	typename LBM_TYPE::T_TRAITS::idx x,
	typename LBM_TYPE::T_TRAITS::idx y,
	typename LBM_TYPE::T_TRAITS::idx z
)
{
	using map_t = typename LBM_TYPE::T_TRAITS::map_t;
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	using dreal = typename LBM_TYPE::T_TRAITS::dreal;

	idx gi = POS(x, y, z, lbm.X, lbm.Y, lbm.Z);
	map_t gi_map = lbm.map(gi);

	KernelStruct<dreal> KS;
//	for (int i=0;i<27;i++) KS.f[i] = lbm.lat[Fxyz(i,x,y,z,lbm.X,lbm.Y,lbm.Z)];
	for (int i=0;i<27;i++) KS.f[i] = lbm.hfs[df_cur][Fxyz((idx)i,x,y,z,lbm.X,lbm.Y,lbm.Z)];

	// copy quantities
	KS.lbmViscosity = lbm.data.lbmViscosity;
	KS.fx=lbm.data.fx;
	KS.fy=lbm.data.fy;
	KS.fz=lbm.data.fz;

	// compute Density & Velocity
	if (LBM_BC::isComputeDensityAndVelocity(gi_map))
		LBM_TYPE::computeDensityAndVelocity(KS);

	MACRO::outputMacro(lbm.hmacro,lbm.X*lbm.Y*lbm.Z,KS,gi);
}


//template<typename L, typename M, typename LBM_DATA>
template <
	typename LBM_TYPE,
	typename STREAMING,
	typename MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
#ifdef USE_CUDA
__global__ void cudaLBMComputeVelocitiesStar(LBM_DATA SD)
#else
void LBMComputeVelocitiesStar(
	LBM_DATA SD,
	typename LBM_TYPE::T_TRAITS::idx x,
	typename LBM_TYPE::T_TRAITS::idx y,
	typename LBM_TYPE::T_TRAITS::idx z
)
#endif
{
	using dreal = typename LBM_TYPE::T_TRAITS::dreal;
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	using map_t = typename LBM_TYPE::T_TRAITS::map_t;

	#ifdef USE_CUDA
	idx x = threadIdx.x + blockIdx.x * blockDim.x;
	idx y = threadIdx.y + blockIdx.y * blockDim.y;
	idx z = threadIdx.z + blockIdx.z * blockDim.z;
	#endif
	idx gi = POS(x, y, z, SD.X, SD.Y, SD.Z);
	map_t gi_map = SD.dmap[gi];

	KernelStruct<dreal> KS;

	// copy quantities
	MACRO::copyQuantities(SD,KS,gi);

	idx xp,xm,yp,ym,zp,zm;
	if (LBM_BC::isPeriodic(gi_map))
	{
		// handle overlaps between GPUs
		xp = (!SD.overlap_right && x == SD.X-1) ? 0 : (x+1);
		xm = (!SD.overlap_left && x == 0) ? (SD.X-1) : (x-1);
		yp = (y == SD.Y-1) ? 0 : (y+1);
		ym = (y == 0) ? (SD.Y-1) : (y-1);
		zp = (z == SD.Z-1) ? 0 : (z+1);
		zm = (z == 0) ? (SD.Z-1) : (z-1);
	} else {
		// handle overlaps between GPUs
		xp = (SD.overlap_right) ? x+1 : MIN(x+1, SD.X-1);
		xm = (SD.overlap_left) ? x-1 : MAX(x-1,0);
		yp = MIN(y+1, SD.Y-1);
		ym = MAX(y-1,0);
		zp = MIN(z+1, SD.Z-1);
		zm = MAX(z-1,0);
	}

	// Streaming
	if (LBM_BC::isStreaming(gi_map))
		STREAMING::streaming(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);

	KS.fx=0;
	KS.fy=0;
	KS.fz=0;

	// compute Density & Velocity
	if (LBM_BC::isComputeDensityAndVelocity(gi_map))
		LBM_TYPE::computeDensityAndVelocity(KS);


	MACRO::outputMacro(SD.dmacro,SD.X*SD.Y*SD.Z,KS,gi);
}

//template<typename L, typename M, typename LBM_DATA>
template <
	typename LBM_TYPE,
	typename STREAMING,
	typename MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
#ifdef USE_CUDA
__global__ void cudaLBMComputeVelocitiesStarAndZeroForce(LBM_DATA SD)
#else
void LBMComputeVelocitiesStarAndZeroForce(
	LBM_DATA SD,
	typename LBM_TYPE::T_TRAITS::idx x,
	typename LBM_TYPE::T_TRAITS::idx y,
	typename LBM_TYPE::T_TRAITS::idx z
)
#endif
{
	using dreal = typename LBM_TYPE::T_TRAITS::dreal;
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	using map_t = typename LBM_TYPE::T_TRAITS::map_t;

 	#ifdef USE_CUDA
	idx x = threadIdx.x + blockIdx.x * blockDim.x;
	idx y = threadIdx.y + blockIdx.y * blockDim.y;
	idx z = threadIdx.z + blockIdx.z * blockDim.z;
	#endif
	idx gi = POS(x, y, z, SD.X, SD.Y, SD.Z);
	map_t gi_map = SD.dmap[gi];

	KernelStruct<dreal> KS;

	// copy quantities
	MACRO::copyQuantities(SD,KS,gi);

	idx xp,xm,yp,ym,zp,zm;
	if (LBM_BC::isPeriodic(gi_map))
	{
		// handle overlaps between GPUs
		xp = (!SD.overlap_right && x == SD.X-1) ? 0 : (x+1);
		xm = (!SD.overlap_left && x == 0) ? (SD.X-1) : (x-1);
		yp = (y == SD.Y-1) ? 0 : (y+1);
		ym = (y == 0) ? (SD.Y-1) : (y-1);
		zp = (z == SD.Z-1) ? 0 : (z+1);
		zm = (z == 0) ? (SD.Z-1) : (z-1);
	} else {
		// handle overlaps between GPUs
		xp = (SD.overlap_right) ? x+1 : MIN(x+1, SD.X-1);
		xm = (SD.overlap_left) ? x-1 : MAX(x-1,0);
		yp = MIN(y+1, SD.Y-1);
		ym = MAX(y-1,0);
		zp = MIN(z+1, SD.Z-1);
		zm = MAX(z-1,0);
	}

	// Streaming
	if (LBM_BC::isStreaming(gi_map))
		STREAMING::streaming(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);

	KS.fx=0;
	KS.fy=0;
	KS.fz=0;

	// compute Density & Velocity
	if (LBM_BC::isComputeDensityAndVelocity(gi_map))
		LBM_TYPE::computeDensityAndVelocity(KS);

	MACRO::outputMacro(SD.dmacro,SD.X*SD.Y*SD.Z,KS,gi);
	// reset forces
	MACRO::zeroForces(SD.dmacro,SD.X*SD.Y*SD.Z,gi);
}


template <
	typename STATE,
	typename LBM
>
void SimUpdate(STATE& state, LBM& lbm)
{
	using MACRO = typename STATE::T_MACRO;
	using CPU_MACRO = typename STATE::T_CPU_MACRO;
	using LBM_TYPE = typename STATE::T_LBM_TYPE;
	using TRAITS = typename LBM_TYPE::T_TRAITS;
	using LBM_DATA = typename STATE::T_LBM_DATA;
	using LBM_BC = typename STATE::T_LBM_BC;

	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	using STREAMING = LBM_STREAMING< TRAITS >;

	// debug
	if (lbm.data.lbmViscosity == 0) {
		state.log("error: LBM viscosity is 0");
		state.lbm.terminate = true;
		return;
	}
	
	
	#ifdef USE_CUDA
		checkCudaDevice;
		dim3 blockSize(1, lbm.block_size, 1);
		dim3 gridSize(lbm.local_X, lbm.Y/lbm.block_size, lbm.Z);

		// check for PEBKAC problem existing between keyboard and chair
		if (gridSize.y * lbm.block_size != lbm.Y) {
			state.log("error: lbm.Y (which is %d) is not aligned to a multiple of the block size (which is %d)", lbm.Y, lbm.block_size);
			state.lbm.terminate = true;
			return;
		}
	#endif
	
	// flags
	bool doComputeVelocitiesStar=false;
	bool doCopyQuantitiesStarToHost=false;
	bool doZeroForceOnDevice=false;
	bool doZeroForceOnHost=false;
	bool doComputeLagrangePhysics=false;
	bool doCopyForceToDevice=false;
    
	// determine global flags
	// NOTE: all Lagrangian points are assumed to be on the first GPU
	if (lbm.offset_X == 0 && state.FF.size() > 0) 
	{
		doComputeLagrangePhysics=true;
		for (int i=0;i<state.FF.size();i++)
		if (state.FF[i].implicitWuShuForcing)
		{
			doComputeVelocitiesStar=true;
			switch (state.FF[i].ws_compute)
			{
				case ws_computeCPU:
					doCopyQuantitiesStarToHost=true;
					doZeroForceOnHost=true;
					doCopyForceToDevice=true;
				case ws_computeCPU_TNL:
					doCopyQuantitiesStarToHost=true;
					doCopyForceToDevice=true;
					break;
				case ws_computeGPU_TNL:
				case ws_computeHybrid_TNL:
				case ws_computeGPU_CUSPARSE:
				case ws_computeHybrid_CUSPARSE:
					doZeroForceOnDevice=true;
					break;
			}
		}
	}


	if (doComputeVelocitiesStar)
	{
		#ifdef USE_CUDA
			if (doZeroForceOnDevice)
				cudaLBMComputeVelocitiesStarAndZeroForce< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC ><<<gridSize, blockSize>>>(lbm.data);
			else
				cudaLBMComputeVelocitiesStar< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC ><<<gridSize, blockSize>>>(lbm.data);
			checkCudaDevice;
		#else
			#pragma omp parallel for schedule(static) collapse(2)
			for (idx x=0;x<lbm.X;x++)
			for (idx z=0;z<lbm.Z;z++)
			for (idx y=0;y<lbm.Y;y++)
			if (doZeroForceOnDevice)
				LBMComputeVelocitiesStarAndZeroForce< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC >(lbm.data, x, y, z);
			else
				LBMComputeVelocitiesStar< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC >(lbm.data, x, y, z);
		#endif
		if (doCopyQuantitiesStarToHost)
		{
			lbm.copyMacroToHost();
		}
	}


	// reset lattice force vectors dfx and dfy
	if (doZeroForceOnHost)
	{
		lbm.resetForces();
	}

//	state.log("core.h state.computeAllLagrangeForces() start");
	if (doComputeLagrangePhysics)
	{
		state.computeAllLagrangeForces();
	}
//	state.log("core.h state.computeAllLagrangeForces() done");

	if (doCopyForceToDevice)
	{
		lbm.copyForcesToDevice();
	}


	#ifdef USE_CUDA
		if (MACRO::use_kernelWorker) cudaMacroWorker< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC><<<gridSize, blockSize>>>(lbm.data);
		cudaLBMKernel< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC><<<gridSize, blockSize>>>(lbm.data);
//		if (MACRO::use_kernelWorker)
		checkCudaDevice;
	#else
		#pragma omp parallel for schedule(static) collapse(2)
		for (idx GX=0;GX<lbm.X;GX++)
		for (idx GZ=0;GZ<lbm.Z;GZ++)
		for (idx GY=0;GY<lbm.Y;GY++)
		{
			if (MACRO::use_kernelWorker) MacroWorker< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC>(lbm.data, GX, GY, GZ);
			LBMKernel< LBM_TYPE, STREAMING, MACRO, LBM_DATA, LBM_BC>(lbm.data, GX, GY, GZ);
		}
	#endif
	
	lbm.iterations++;

	bool doit=false; 
	for (int c=0;c<MAX_COUNTER;c++) if (c!=PRINT && c!=SAVESTATE) if (state.cnt[c].action(lbm.physTime())) doit = true;
	if (doit)
	{
		lbm.copyMacroToHost();
		if (CPU_MACRO::N>0) lbm.copyDFsToHost(lbm.dfs[df_out], lbm.hfs[df_out]); // to be able to compute rho, vx, vy, vz etc... based on DFs on CPU to save GPU memory FIXME may not work with ESOTWIST
		#ifdef USE_CUDA
		checkCudaDevice;
		#endif
	}
}

// called just once after the SimUpdate - only global stuff should be here
template <
	typename STATE
>
void AfterSimUpdate(STATE& state, timespec& t1, timespec& t2, int& lbmPrevIterations)
{
	typename STATE::T_LBM& lbm = state.lbm;

	#ifdef USE_CUDA
	// synchronization is not necessary for correctness, only to get correct MLUPS
	bool doit=false; 
	for (int c=0;c<MAX_COUNTER;c++) if (state.cnt[c].action(lbm.physTime())) doit = true;
	if (doit)
	{
		cudaDeviceSynchronize();
		checkCudaDevice;
	}
	#endif

	bool write_info=false;

	if (state.cnt[VTK1D].action(lbm.physTime()) || 
	    state.cnt[VTK2D].action(lbm.physTime()) ||
	    state.cnt[VTK3D].action(lbm.physTime()) ||
	    state.cnt[VTK3DCUT].action(lbm.physTime()) ||
	    state.cnt[PROBE1].action(lbm.physTime()) ||
	    state.cnt[PROBE2].action(lbm.physTime()) ||
	    state.cnt[PROBE3].action(lbm.physTime())
	    )
	{
		// common copy
		state.lbm.computeCPUMacroFromLat();
		// probe1
		if (state.cnt[PROBE1].action(lbm.physTime()))
		{
			state.probe1();
			state.cnt[PROBE1].count++;
		}
		// probe2
		if (state.cnt[PROBE2].action(lbm.physTime()))
		{
			state.probe2();
			state.cnt[PROBE2].count++;
		}
		// probe3
		if (state.cnt[PROBE3].action(lbm.physTime()))
		{
			state.probe3();
			state.cnt[PROBE3].count++;
		}
		// 3D VTK
		if (state.cnt[VTK3D].action(lbm.physTime()))
		{
			state.writeVTKs_3D();
			state.cnt[VTK3D].count++;
		}
		// 3D VTK CUT
		if (state.cnt[VTK3DCUT].action(lbm.physTime()))
		{
			state.writeVTKs_3Dcut();
			state.cnt[VTK3DCUT].count++;
		}
		// 2D VTK
		if (state.cnt[VTK2D].action(lbm.physTime()))
		{
			state.writeVTKs_2D();
			state.cnt[VTK2D].count++;
		}
		// 1D VTK
		if (state.cnt[VTK1D].action(lbm.physTime()))
		{
			state.writeVTKs_1D();
			state.cnt[VTK1D].count++;
		}
		write_info = true;
	}

	if (state.cnt[PRINT].action(lbm.physTime()))
	{
		write_info = true;
		state.cnt[PRINT].count++;
	}

/*
	if (state.cnt[SAVESTATE].action(lbm.physTime()))
	{
        printf("provadim saveState\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
		write_info = true;
		state.saveState();
		state.cnt[SAVESTATE].count++;
	}
*/	


	// statReset is called after all probes and VTK output
//	if (state.statResetPeriod > 0 && lbm.physTime() >= state.statResetCount * state.statResetPeriod)
	if (state.cnt[STAT_RESET].action(lbm.physTime())) state.statReset();
	if (state.cnt[STAT2_RESET].action(lbm.physTime())) state.stat2Reset();

	clock_gettime(CLOCK_REALTIME, &t2);

	if (lbm.iterations > 1)
//	if (write_info || (state.printIter > 0 && lbm.iterations % state.printIter == 0) )
	if (write_info)
	{
		long timediff = (t2.tv_sec - t1.tv_sec) * 1000000000 + (t2.tv_nsec - t1.tv_nsec);
		// to avoid numerical errors - split LUPS computation in two parts
		double LUPS = lbm.iterations - lbmPrevIterations;
		LUPS *= lbm.X * lbm.Y * lbm.Z * 1000000000.0 / timediff;
		write_info = true;
		clock_gettime(CLOCK_REALTIME, &t1);
		lbmPrevIterations=lbm.iterations;

		if (state.verbosity>0)
		{
			state.log("MLUPS=%.1f iter=%d t=%1.3fs dt=%1.2e lbmVisc=%1.2e",
				LUPS * 1e-6,
				lbm.iterations,
				lbm.physTime(),
				lbm.physDt,
				lbm.lbmViscosity()
			);
		}
	}
}

// called multiple times after AfterSimUpdate for each GPU - local stuff should be handled here
template <
	typename STATE, 
	typename LBM
>
void AfterSimUpdate_local(STATE& state, LBM& lbm, int cpu_thread_id = 0)
{
	// copy macro from host to device after reset
	if (state.cnt[STAT_RESET].action(lbm.physTime()))
	{
		lbm.copyMacroToDevice();
		if (cpu_thread_id == 0)
			state.cnt[STAT_RESET].count++;
	}
	if (state.cnt[STAT2_RESET].action(lbm.physTime()))
	{
		lbm.copyMacroToDevice();
		if (cpu_thread_id == 0)
			state.cnt[STAT2_RESET].count++;
	}
	
	//check wall time
	if (state.cnt[SAVESTATE].action(lbm.physTime()) || state.wallTimeReached())
	{
		lbm.copyMacroToHost();
		lbm.copyDFsToHost();
		lbm.copyMapToHost();
		// TODO: copy all quantities to CPU from lbm
//		lbm.copyMacroToDevice();
//		if (cpu_thread_id == 0)
//			state.cnt[STAT2_RESET].count++;
	}

}

template <
	typename STATE
>
bool estimateMemoryDemands(STATE& state)
{
	using MACRO = typename STATE::T_MACRO;
	using LBM_TYPE = typename STATE::T_LBM_TYPE;
	using LBM_DATA = typename STATE::T_LBM_DATA;
	using LBM = typename STATE::T_LBM;
	using TRAITS = typename STATE::T_TRAITS;
	using LBM_BC = typename STATE::T_LBM_BC;
	
	using idx = typename TRAITS::idx;
	using idx = typename TRAITS::idx;
	using real = typename TRAITS::real;
	using dreal = typename TRAITS::dreal;
	using map_t = typename TRAITS::map_t;

	typename STATE::T_LBM& lbm = state.lbm; // only a reference to state.lbm is stored to lbm

	long long memDFs = lbm.X*lbm.Y*lbm.Z*27*sizeof(dreal);
	long long memMacro = lbm.X*lbm.Y*lbm.Z*sizeof(dreal)*MACRO::N;
	long long memMap = lbm.X*lbm.Y*lbm.Z*sizeof(map_t);
	long long CPUavail = sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE);
	long long GPUavail = 0;
	long long GPUtotal = 0;
	long long GPUtotal_hw = 0;
	long long CPUtotal = memMacro + memMap + DFMAX*memDFs; 
	long long CPUDFs = DFMAX*memDFs;
	#ifdef USE_CUDA
	GPUavail = 0;
	GPUtotal_hw =0;
	GPUtotal += DFMAX*memDFs + memMacro + memMap;
//	CPUDFs = 0;

	// get number of CUDA GPUs
	int num_gpus=0;
	cudaGetDeviceCount(&num_gpus);

	// display CPU and GPU configuration
	state.log("number of CUDA devices:\t%d", num_gpus);
	for (int i = 0; i < num_gpus; i++)
	{
		cudaDeviceProp dprop;
		cudaGetDeviceProperties(&dprop, i);
		state.log("   %d: %s", i, dprop.name);
		cudaSetDevice(i);
		size_t free=0, total=0;
		cudaMemGetInfo(&free, &total);
		GPUavail += free;
		GPUtotal_hw += total;
	}

	#else
//	CPUtotal += CPUDFs;
	#endif
	
	state.log("Memory budget analysis / estimation");
	state.log("CPU RAM for DFs:   %ld MiB", (long)(CPUDFs/1024/1024));
//	state.log("CPU RAM for lat:   %ld MiB", (long)(memDFs/1024/1024));
	state.log("CPU RAM for map:   %ld MiB", (long)(memMap/1024/1024));
	state.log("CPU RAM for macro: %ld MiB", (long)(memMacro/1024/1024));
	state.log("TOTAL CPU RAM %ld MiB estimated needed, %ld MiB available (%6.4f%%)", (long)(CPUtotal/1024/1024), (long)(CPUavail/1024/1024), (double)(100.0*CPUtotal/CPUavail));
	#ifdef USE_CUDA
	state.log("GPU RAM for DFs:   %ld MiB", (long)(2*memDFs/1024/1024));
	state.log("GPU RAM for map:   %ld MiB", (long)(memMap/1024/1024));
	state.log("GPU RAM for macro: %ld MiB", (long)(memMacro/1024/1024));
	state.log("TOTAL GPU RAM %ld MiB estimated needed, %ld MiB available (%6.4f%%), total GPU RAM: %ld MiB", (long)(GPUtotal/1024/1024), (long)(GPUavail/1024/1024), (double)(100.0*GPUtotal/GPUavail), (long)(GPUtotal_hw/1024/1024));
	if (GPUavail <= GPUtotal) return false;
	#endif
	if (CPUavail <= CPUtotal) return false;
	return true;
}


template <
	typename STATE
>
void execute(STATE& state)
{
	using MACRO = typename STATE::T_MACRO;
	using LBM_TYPE = typename STATE::T_LBM_TYPE;
	using LBM_DATA = typename STATE::T_LBM_DATA;
	using LBM = typename STATE::T_LBM;
	using TRAITS = typename STATE::T_TRAITS;
	using LBM_BC = typename STATE::T_LBM_BC;
	
	using idx = typename TRAITS::idx;
	using idx = typename TRAITS::idx;
	
	// for parallel vtk write
	omp_set_nested(1);

	if (!estimateMemoryDemands(state))
	{
		state.log("Not enough memory available (CPU or GPU). [disable this check in lbm3d/core.h -> execute]");
		return;
	}

	typename STATE::T_LBM& lbm = state.lbm; // only a reference to state.lbm is stored to lbm
	
	
	// check for loadState
	// reset counters
	for (int c=0;c<MAX_COUNTER;c++) state.cnt[c].count=0;
	lbm.iterations=0;
	int lbmPrevIterations=0;

	state.reset(); // initial condition --> lat or loadState
    

	#ifdef USE_CUDA
	// create CUDA streams
	for (int i = 0; i < max_cuda_streams; i++) {
		cudaStreamCreate(&cuda_streams[i]);
	}
	checkCudaDevice;
	#endif

	// initialize macroscopic quantities 
	#pragma omp parallel for schedule(static) collapse(2)
	for (idx GX=0;GX<lbm.X;GX++)
	for (idx GZ=0;GZ<lbm.Z;GZ++)
	for (idx GY=0;GY<lbm.Y;GY++)
		LBMKernelInit< LBM_TYPE, LBM_BC, MACRO, LBM>(lbm, GX, GY, GZ);

	#ifdef USE_CUDA
	checkCudaDevice;
	#endif

	// must be volatile, otherwise OpenMP deadlocks at the end of the main loop
	volatile bool quit = false;

	struct timespec t1, t2;
	clock_gettime(CLOCK_REALTIME, &t1);

	state.log("\nSTART: simulation LBM:%s dimensions %d x %d x %d lbmVisc %e physDl %e physDt %e", STATE::T_LBM_TYPE::id, lbm.X, lbm.Y, lbm.Z, lbm.lbmViscosity(), lbm.physDl, lbm.physDt);

	#ifdef USE_CUDA
	int num_gpus = 0;
	if (lbm.use_multiple_gpus)
	{
		// get number of CUDA GPUs
		cudaGetDeviceCount(&num_gpus);

		if (num_gpus < 1)
		{
			state.log("ERROR: no CUDA capable devices were detected");
			return;
		}

		// display CPU and GPU configuration
		state.log("number of CPU cores:\t%d", omp_get_num_procs());
		state.log("number of CUDA devices:\t%d", num_gpus);
		for (int i = 0; i < num_gpus; i++)
		{
			cudaDeviceProp dprop;
			cudaGetDeviceProperties(&dprop, i);
			state.log("   %d: %s", i, dprop.name);
		}
	
		// hot fix FIXME
		if (num_gpus < 2)
		{
			state.log("warning: multi GPU on a single GPU, num_gpus %d", num_gpus);
//			lbm.use_multiple_gpus = false;
		}
	}
	#else
	lbm.use_multiple_gpus = false;
	#endif



	#ifdef USE_CUDA
	if (lbm.use_multiple_gpus)
	{
		////////////////////////////////////////////////////////////////
		// run as many CPU threads as there are CUDA devices
		//   each CPU thread controls a different device, processing its
		//   portion of the data.  It's possible to use more CPU threads
		//   than there are CUDA devices, in which case several CPU
		//   threads will be allocating resources and launching kernels
		//   on the same device.  For example, try omp_set_num_threads(2*num_gpus);
		//
		int num_cpu_threads = num_gpus;  // create as many CPU threads as there are CUDA devices
//		int num_cpu_threads = 2*num_gpus;  // create twice as many CPU threads as there are CUDA devices

		// array for sharing LBM instances for the synchronization of overlaps
		typename STATE::T_LBM* lbms[num_cpu_threads];
		int gpu_ids[num_cpu_threads];

		omp_set_num_threads(num_cpu_threads);
		#pragma omp parallel
		{
			int cpu_thread_id = omp_get_thread_num();
			int gpu_id = cpu_thread_id % num_gpus;   // "% num_gpus" allows more CPU threads than GPU devices
			gpu_ids[cpu_thread_id] = gpu_id;

			#pragma omp barrier

			// neighbors
			int cpu_id_left = (cpu_thread_id - 1 + num_cpu_threads) % num_cpu_threads;
			int gpu_id_left = gpu_ids[cpu_id_left];
			int cpu_id_right = (cpu_thread_id + 1) % num_cpu_threads;
			int gpu_id_right = gpu_ids[cpu_id_right];

			// set and check the CUDA device for this CPU thread
			cudaSetDevice(gpu_id);
			checkCudaDevice;

			// local X size
			idx local_X = lbm.X / num_cpu_threads;
			idx offset_X = local_X * cpu_thread_id;
			if (cpu_thread_id == num_cpu_threads - 1)
				local_X += lbm.X % num_cpu_threads;

			// create a local lbm instance
			typename STATE::T_LBM local_lbm = lbm;
			lbms[cpu_thread_id] = &local_lbm;
			
			// maybe useless/not needed: TODO check
			local_lbm.hmap = lbm.hmap;
			for (int dfty=0;dfty<DFMAX;dfty++) local_lbm.hfs[dfty]=lbm.hfs[dfty];
//			local_lbm.lat = lbm.lat;
			local_lbm.wall = lbm.wall;
			local_lbm.hmacro = lbm.hmacro;
			local_lbm.cpumacro = lbm.cpumacro;

			// set offset and overlaps
			local_lbm.local_X = local_X;
			local_lbm.offset_X = offset_X;
			// NOTE: extra copy can be avoided iff the x-normal boundaries are not periodic
//			if (cpu_thread_id > 0)
			local_lbm.overlap_left = 1;
			// NOTE: extra copy can be avoided iff the x-normal boundaries are not periodic
//			if (cpu_thread_id < num_cpu_threads - 1)
			local_lbm.overlap_right = 1;

			#pragma omp barrier

			// setup peer access
			if (local_lbm.overlap_left) {
				int access;
				cudaDeviceCanAccessPeer(&access, gpu_id, gpu_id_left);
				state.log("Peer access from GPU %d to GPU %d: %s", gpu_id, gpu_id_left, access ? "possible" : "NOT possible");
				if (access)
					cudaDeviceEnablePeerAccess(gpu_id_left, 0);
			}
			if (local_lbm.overlap_right && gpu_id_left != gpu_id_right) {
				int access;
				cudaDeviceCanAccessPeer(&access, gpu_id, gpu_id_right);
				state.log("Peer access from GPU %d to GPU %d: %s", gpu_id, gpu_id_right, access ? "possible" : "NOT possible");
				if (access)
					cudaDeviceEnablePeerAccess(gpu_id_right, 0);
			}
			checkCudaDevice;

			local_lbm.allocateDeviceData();
			local_lbm.copyMapToDevice();
			local_lbm.copyDFsToDevice();
			local_lbm.copyMacroToDevice();  // important when a state has been loaded
			checkCudaDevice;
			#pragma omp barrier

			if (cpu_thread_id == 0)
			{
				// make snapshot for the initial condition
				AfterSimUpdate(state, t1, t2, lbmPrevIterations);
			}
			#pragma omp barrier
			AfterSimUpdate_local(state, local_lbm, cpu_thread_id);

			// sync Macro (for Initial condition only)
			// synchronize macro
			if (MACRO::use_syncMacro)
			{
				if (local_lbm.overlap_right)
					local_lbm.synchronizeMacroToRight(gpu_id, gpu_id_right, *lbms[cpu_id_right]);
				if (local_lbm.overlap_left)
					local_lbm.synchronizeMacroToLeft(gpu_id, gpu_id_left, *lbms[cpu_id_left]);
			}

			#pragma omp barrier

			while (!quit)
			{
				// synchronize overlaps
				if (local_lbm.overlap_right)
					local_lbm.synchronizeDFsToRight(gpu_id, gpu_id_right, *lbms[cpu_id_right]);
				if (local_lbm.overlap_left)
					local_lbm.synchronizeDFsToLeft(gpu_id, gpu_id_left, *lbms[cpu_id_left]);

				cudaDeviceSynchronize();
				#pragma omp barrier

				// update kernel data (viscosity, swap df1 and df2)
				local_lbm.updateKernelData();
				state.updateKernelVelocities(local_lbm);

				// run LBM kernel
				SimUpdate(state, local_lbm);

				cudaDeviceSynchronize();
				#pragma omp barrier

				// synchronize macro		// FIXME: order
				if (MACRO::use_syncMacro)
				{
					if (local_lbm.overlap_right)
						local_lbm.synchronizeMacroToRight(gpu_id, gpu_id_right, *lbms[cpu_id_right]);
					if (local_lbm.overlap_left)
						local_lbm.synchronizeMacroToLeft(gpu_id, gpu_id_left, *lbms[cpu_id_left]);
				}
				#pragma omp barrier

				if (cpu_thread_id == 0)
				{
					// sync iterations counter
					lbm.iterations = local_lbm.iterations;
					// post-processing: snapshots etc.
					AfterSimUpdate(state, t1, t2, lbmPrevIterations);
				}
				#pragma omp barrier
				AfterSimUpdate_local(state, local_lbm, cpu_thread_id);

				cudaDeviceSynchronize();
				#pragma omp barrier
				// only the master thread handles time stepping
				if (cpu_thread_id == 0)
				{
					//check wall time
					if (state.wallTimeReached())
					{
						// copy data to CPU (if needed)
						state.saveState(true);
						quit = true;
					}
					
					// check savestate
					if (state.cnt[SAVESTATE].action(lbm.physTime()))
					{
						state.saveState();
						state.cnt[SAVESTATE].count++;
					}

					if (lbm.physTime() > lbm.physFinalTime) quit = true;
					if (lbm.quit())
					{
						state.log("terminate flag triggered");
						quit=true;
					}
				}

				// avoid deadlock when one thread sets quit=true
				#pragma omp barrier
			}

			// disable peer access to make CUDA happy (otherwise the next process might get "peer access is already enabled" error)
			if (local_lbm.overlap_left) {
				int access;
				cudaDeviceCanAccessPeer(&access, gpu_id, gpu_id_left);
				if (access)
					cudaDeviceDisablePeerAccess(gpu_id_left);
			}
			if (local_lbm.overlap_right && gpu_id_left != gpu_id_right) {
				int access;
				cudaDeviceCanAccessPeer(&access, gpu_id, gpu_id_right);
				if (access)
					cudaDeviceDisablePeerAccess(gpu_id_right);
			}
			checkCudaDevice;

//			cudaProfilerStop();

			// avoid double free
			local_lbm.hmap = NULL;
			//local_lbm.lat = NULL;
			for (int dfty=0;dfty<DFMAX;dfty++) local_lbm.hfs[dfty]=NULL;
			local_lbm.wall = NULL;
			local_lbm.hmacro = NULL;
			local_lbm.cpumacro = NULL;
		}
	}
	else
	{
		lbm.allocateDeviceData();
		lbm.copyMapToDevice();
		lbm.copyDFsToDevice();
		lbm.copyMacroToDevice();  // important when a state has been loaded
		checkCudaDevice;
		// make snapshot for the initial condition
		AfterSimUpdate(state, t1, t2, lbmPrevIterations);
		AfterSimUpdate_local(state, lbm);
		while (!quit)
		{
			state.lbm.updateKernelData();
			state.updateKernelVelocities(state.lbm);
			SimUpdate(state, lbm);
			AfterSimUpdate(state, t1, t2, lbmPrevIterations);
			AfterSimUpdate_local(state, lbm);

			//check wall time
			if (state.wallTimeReached())
			{
				// copy data to CPU (if needed)
				state.saveState(true);
				quit = true;
			}
			if (state.cnt[SAVESTATE].action(lbm.physTime()))
			{
				// save state -- all quantities need to be copied to CPU
				state.saveState();
				state.cnt[SAVESTATE].count++;
			}
            
			if (lbm.physTime() > lbm.physFinalTime) quit = true;
			if (lbm.quit())
			{
				state.log("terminate flag triggered");
				quit=true;
			}
		}
	}
	#else
	lbm.allocateDeviceData();
	lbm.copyMapToDevice();  // not needed here
	lbm.copyDFsToDevice();  // not needed here
	lbm.copyMacroToDevice(); // not needed here

	// make snapshot for the initial condition
	AfterSimUpdate(state, t1, t2, lbmPrevIterations);
	AfterSimUpdate_local(state, lbm);
	
	while (!quit)
	{
		state.lbm.updateKernelData();
		state.updateKernelVelocities(state.lbm);
		SimUpdate(state, lbm);
		AfterSimUpdate(state, t1, t2, lbmPrevIterations);
		AfterSimUpdate_local(state, lbm);

		if (state.wallTimeReached())
		{
			state.saveState(true);
			quit = true;
		}

		if (state.cnt[SAVESTATE].action(lbm.physTime()))
		{
			// save state -- all quantities need to be copied to CPU
			state.saveState();
			state.cnt[SAVESTATE].count++;
		}

		if (lbm.physTime() > lbm.physFinalTime) quit = true;
		if (lbm.quit())
		{
			state.log("terminate flag triggered");
			quit=true;
		}
	}
	#endif // USE_CUDA

	#ifdef USE_CUDA
	// destroy CUDA streams
	for (int i = 0; i < max_cuda_streams; i++)
		cudaStreamDestroy(cuda_streams[i]);
	checkCudaDevice;
	#endif
}

