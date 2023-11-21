#ifndef _LBM_EQ_ENTROPIC_H_
#define _LBM_EQ_ENTROPIC_H_

// entropic equilibrium (for KBC) by Robert Straka
//Entropic equilibrium -- somebody should check *FIXME* if this type of equilibrium gives correct 3rd order moment e.g. M^eq_111=rho*u_x*u_y*u_z as
// this is not the case for standard quadratic equilibria! check papers by P.J. Dellar for further details how to handle this cubic defect.

#define _efeq_mmm(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_mmz(vx,vy,vz) (n1o6*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_mmp(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_mzm(vx,vy,vz) (n1o6*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_mzz(vx,vy,vz) (n1o6*n2o3*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx)))
#define _efeq_mzp(vx,vy,vz) (n1o6*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_mpm(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_mpz(vx,vy,vz) (n1o6*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_mpp(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zmm(vx,vy,vz) (n2o3*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zmz(vx,vy,vz) (n2o3*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_zmp(vx,vy,vz) (n2o3*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zzm(vx,vy,vz) (n2o3*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zzz(vx,vy,vz) (n2o3*n2o3*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz)))
#define _efeq_zzp(vx,vy,vz) (n2o3*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zpm(vx,vy,vz) (n2o3*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_zpz(vx,vy,vz) (n2o3*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_zpp(vx,vy,vz) (n2o3*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_pmm(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_pmz(vx,vy,vz) (n1o6*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_pmp(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_pzm(vx,vy,vz) (n1o6*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_pzz(vx,vy,vz) (n1o6*n2o3*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx)))
#define _efeq_pzp(vx,vy,vz) (n1o6*n2o3*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_ppm(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*no1/((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))
#define _efeq_ppz(vx,vy,vz) (n1o6*n1o6*n2o3*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy)))
#define _efeq_ppp(vx,vy,vz) (n1o6*n1o6*n1o6*(no2 - sqrt(no1 + no3*vx*vx))*(no2 - sqrt(no1 + no3*vy*vy))*(no2 - sqrt(no1 + no3*vz*vz))*((no2*vx + sqrt(no1 + no3*vx*vx))/(no1 - vx))*((no2*vy + sqrt(no1 + no3*vy*vy))/(no1 - vy))*((no2*vz + sqrt(no1 + no3*vz*vz))/(no1 - vz)))


template < typename TRAITS >
struct LBM_EQ_ENTROPIC
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	CUDA_HOSTDEV static dreal feq_zzz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zzz(vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_pzz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pzz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mzz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zpz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zmz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zzp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zzp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zzm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zzm(vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_ppz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_ppz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pmz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mpz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmz(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mmz(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pzp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pzp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mzm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pzm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pzm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mzp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mzp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zpp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zpm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zpm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zmp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_zmm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_zmm(vx, vy, vz);}

	CUDA_HOSTDEV static dreal feq_ppp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_ppp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mmm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_ppm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_ppm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pmp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mpp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mpm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mpm(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_mmp(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_mmp(vx, vy, vz);}
	CUDA_HOSTDEV static dreal feq_pmm(dreal rho, dreal vx, dreal vy, dreal vz) {return rho*_efeq_pmm(vx, vy, vz);}
};

#endif

