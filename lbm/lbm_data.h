#pragma once

template < typename TRAITS >
struct LBM_Data
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;
	using map_t = typename TRAITS::map_t;

//	idx X, Y, Z, offset_X, X_plus_overlaps; // FIXME
	idx X, Y, Z;//, offset_X;
	int overlap_left=0, overlap_right=0;
	int macro_overlap_left=0, macro_overlap_right=0;

	// homogeneous force field
	dreal fx = 0;
	dreal fy = 0;
	dreal fz = 0;

	dreal lbmViscosity;
	dreal *dfs[DFMAX];
	dreal *dmacro;
	map_t *dmap;

	CUDA_HOSTDEV dreal& df(int type, int q, idx x, idx y, idx z) { return dfs[type][Fxyz(q,x+overlap_left,y,z,X+overlap_left+overlap_right,Y,Z)]; }
};


template < typename TRAITS >
struct LBM_Data_ConstInflow : LBM_Data<TRAITS>
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	dreal inflow_rho=no1;
	dreal inflow_vx=0;
	dreal inflow_vy=0;
	dreal inflow_vz=0;
	CUDA_HOSTDEV void inflow(KernelStruct<dreal> &KS, idx x, idx y, idx z)
	{
		KS.rho = inflow_rho;
		KS.vx  = inflow_vx;
		KS.vy  = inflow_vy;
		KS.vz  = inflow_vz;
	}
};

template < typename TRAITS >
struct LBM_Data_NoInflow : LBM_Data<TRAITS>
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	CUDA_HOSTDEV void inflow(KernelStruct<dreal> &KS, idx x, idx y, idx z)
	{
		KS.rho = no1;
		KS.vx  = 0;
		KS.vy  = 0;
		KS.vz  = 0;
	}
};
