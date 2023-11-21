#ifndef _LBM_EQ_H_
#define _LBM_EQ_H_

//#define feq(qx, qy, qz, vx, vy, vz) no1 - n3o2 * ((vx)*(vx) + (vy)*(vy) + (vz)*(vz)) + no3*((qx)*(vx) + (qy)*(vy) + (qz)*(vz)) + n9o2*((qx)*(vx) + (qy)*(vy) + (qz)*(vz))*((qx)*(vx) + (qy)*(vy) + (qz)*(vz))

// second order Maxwell-Boltzmann Equilibrium
template < typename TRAITS >
struct LBM_EQ_DEFAULT 
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	CUDA_HOSTDEV_NOINLINE static dreal feq(int qx, int qy, int qz, dreal vx, dreal vy, dreal vz)
	{
		return no1 - n3o2 * (vx*vx + vy*vy + vz*vz) + no3*(qx*vx + qy*vy + qz*vz) + n9o2*(qx*vx + qy*vy + qz*vz)*(qx*vx + qy*vy + qz*vz);
	}


	CUDA_HOSTDEV static dreal feq_zzz(dreal rho, dreal vx, dreal vy, dreal vz) {return n8o27*rho*feq( 0, 0, 0,  vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_pzz(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq( 1, 0, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzz(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq(-1, 0, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpz(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq( 0, 1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmz(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq( 0,-1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zzp(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq( 0, 0, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zzm(dreal rho, dreal vx, dreal vy, dreal vz) {return n2o27*rho*feq( 0, 0,-1, vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_ppz(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 1, 1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmz(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 1,-1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpz(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq(-1, 1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmz(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq(-1,-1, 0, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pzp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 1, 0, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq(-1, 0,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pzm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 1, 0,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq(-1, 0, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 0, 1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 0, 1,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 0,-1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o54*rho*feq( 0,-1,-1, vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_ppp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq( 1, 1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq(-1,-1,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_ppm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq( 1, 1,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq( 1,-1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq(-1, 1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq(-1, 1,-1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmp(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq(-1,-1, 1, vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmm(dreal rho, dreal vx, dreal vy, dreal vz) {return n1o216*rho*feq( 1,-1,-1, vx, vy, vz);}
};

//#undef feq

#endif
