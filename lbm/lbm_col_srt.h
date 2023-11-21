// ]improved BRK (SRT) model by Geier 2017
// for standard DF (no well-conditioned)

#include "lbm_common.h"
template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_SRT : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "SRT"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		const dreal tau = no3*KS.lbmViscosity+n1o2;
		// forcing: vzorce_bgk_force.mw 
		const dreal Smmm = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Szmm = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Spmm = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Smzm = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Szzm = (no3*(-KS.vx*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Spzm = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Smpm = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Szpm = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Sppm = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz))/KS.rho;
		const dreal Smmz = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Szmz = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Spmz = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Smzz = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Szzz = (no3*(-KS.fx*KS.vx-KS.fy*KS.vy-KS.fz*KS.vz))/KS.rho;
		const dreal Spzz = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Smpz = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Szpz = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Sppz = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz))/KS.rho;
		const dreal Smmp = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Szmp = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Spmp = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Smzp = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Szzp = (no3*(-KS.vx*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Spzp = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Smpp = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Szpp = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;
		const dreal Sppp = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz))/KS.rho;

		const dreal locfeq_mmm = LBM_EQ::feq_mmm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mmz = LBM_EQ::feq_mmz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mmp = LBM_EQ::feq_mmp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mzm = LBM_EQ::feq_mzm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mzz = LBM_EQ::feq_mzz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mzp = LBM_EQ::feq_mzp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mpm = LBM_EQ::feq_mpm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mpz = LBM_EQ::feq_mpz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_mpp = LBM_EQ::feq_mpp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zmm = LBM_EQ::feq_zmm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zmz = LBM_EQ::feq_zmz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zmp = LBM_EQ::feq_zmp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zzm = LBM_EQ::feq_zzm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zzz = LBM_EQ::feq_zzz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zzp = LBM_EQ::feq_zzp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zpm = LBM_EQ::feq_zpm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zpz = LBM_EQ::feq_zpz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_zpp = LBM_EQ::feq_zpp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pmm = LBM_EQ::feq_pmm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pmz = LBM_EQ::feq_pmz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pmp = LBM_EQ::feq_pmp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pzm = LBM_EQ::feq_pzm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pzz = LBM_EQ::feq_pzz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_pzp = LBM_EQ::feq_pzp(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_ppm = LBM_EQ::feq_ppm(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_ppz = LBM_EQ::feq_ppz(KS.rho,KS.vx,KS.vy,KS.vz);
		const dreal locfeq_ppp = LBM_EQ::feq_ppp(KS.rho,KS.vx,KS.vy,KS.vz);

		KS.f[mmm] += (locfeq_mmm - KS.f[mmm])/tau + (no1 - n1o2/tau)*Smmm*locfeq_mmm;
		KS.f[mmz] += (locfeq_mmz - KS.f[mmz])/tau + (no1 - n1o2/tau)*Smmz*locfeq_mmz;
		KS.f[mmp] += (locfeq_mmp - KS.f[mmp])/tau + (no1 - n1o2/tau)*Smmp*locfeq_mmp;
		KS.f[mzm] += (locfeq_mzm - KS.f[mzm])/tau + (no1 - n1o2/tau)*Smzm*locfeq_mzm;
		KS.f[mzz] += (locfeq_mzz - KS.f[mzz])/tau + (no1 - n1o2/tau)*Smzz*locfeq_mzz;
		KS.f[mzp] += (locfeq_mzp - KS.f[mzp])/tau + (no1 - n1o2/tau)*Smzp*locfeq_mzp;
		KS.f[mpm] += (locfeq_mpm - KS.f[mpm])/tau + (no1 - n1o2/tau)*Smpm*locfeq_mpm;
		KS.f[mpz] += (locfeq_mpz - KS.f[mpz])/tau + (no1 - n1o2/tau)*Smpz*locfeq_mpz;
		KS.f[mpp] += (locfeq_mpp - KS.f[mpp])/tau + (no1 - n1o2/tau)*Smpp*locfeq_mpp;
		KS.f[zmm] += (locfeq_zmm - KS.f[zmm])/tau + (no1 - n1o2/tau)*Szmm*locfeq_zmm;
		KS.f[zmz] += (locfeq_zmz - KS.f[zmz])/tau + (no1 - n1o2/tau)*Szmz*locfeq_zmz;
		KS.f[zmp] += (locfeq_zmp - KS.f[zmp])/tau + (no1 - n1o2/tau)*Szmp*locfeq_zmp;
		KS.f[zzm] += (locfeq_zzm - KS.f[zzm])/tau + (no1 - n1o2/tau)*Szzm*locfeq_zzm;
		KS.f[zzz] += (locfeq_zzz - KS.f[zzz])/tau + (no1 - n1o2/tau)*Szzz*locfeq_zzz;
		KS.f[zzp] += (locfeq_zzp - KS.f[zzp])/tau + (no1 - n1o2/tau)*Szzp*locfeq_zzp;
		KS.f[zpm] += (locfeq_zpm - KS.f[zpm])/tau + (no1 - n1o2/tau)*Szpm*locfeq_zpm;
		KS.f[zpz] += (locfeq_zpz - KS.f[zpz])/tau + (no1 - n1o2/tau)*Szpz*locfeq_zpz;
		KS.f[zpp] += (locfeq_zpp - KS.f[zpp])/tau + (no1 - n1o2/tau)*Szpp*locfeq_zpp;
		KS.f[pmm] += (locfeq_pmm - KS.f[pmm])/tau + (no1 - n1o2/tau)*Spmm*locfeq_pmm;
		KS.f[pmz] += (locfeq_pmz - KS.f[pmz])/tau + (no1 - n1o2/tau)*Spmz*locfeq_pmz;
		KS.f[pmp] += (locfeq_pmp - KS.f[pmp])/tau + (no1 - n1o2/tau)*Spmp*locfeq_pmp;
		KS.f[pzm] += (locfeq_pzm - KS.f[pzm])/tau + (no1 - n1o2/tau)*Spzm*locfeq_pzm;
		KS.f[pzz] += (locfeq_pzz - KS.f[pzz])/tau + (no1 - n1o2/tau)*Spzz*locfeq_pzz;
		KS.f[pzp] += (locfeq_pzp - KS.f[pzp])/tau + (no1 - n1o2/tau)*Spzp*locfeq_pzp;
		KS.f[ppm] += (locfeq_ppm - KS.f[ppm])/tau + (no1 - n1o2/tau)*Sppm*locfeq_ppm;
		KS.f[ppz] += (locfeq_ppz - KS.f[ppz])/tau + (no1 - n1o2/tau)*Sppz*locfeq_ppz;
		KS.f[ppp] += (locfeq_ppp - KS.f[ppp])/tau + (no1 - n1o2/tau)*Sppp*locfeq_ppp;
	}
};
