// ]improved BRK (SRT) model by Geier 2017
// for standard DF (no well-conditioned)
#include "lbm_common.h"


template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_SRT_MODIF_FORCE : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "SRT_MRT_MODIF_FORCE"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		const dreal tau = no3*KS.lbmViscosity+n1o2;
		// Modif_Force_SRT.mw PE
		const dreal Smmm = (2.0*KS.vx +3.0*KS.vy + 3.0*KS.vz - 1.0)*KS.fx/72. + (3.0*KS.vx +2.0*KS.vy + 3.0*KS.vz - 1.0)*KS.fy/72. + (3.0*KS.vx +3.0*KS.vy + 2.0*KS.vz - 1.)*KS.fz/72.;
		const dreal Szmm = (2.0*KS.vy + 3.0*KS.vz - 1.0)*KS.fy/18. + (3.0*KS.vy + 2.0*KS.vz - 1.0)*KS.fz/18. - KS.fx*KS.vx/18.; 
		const dreal Spmm = (2.0*KS.vx - 3.0*KS.vy - 3.0*KS.vz + 1.0)*KS.fx/72. + (-3.0*KS.vx + 2.0*KS.vy + 3.0*KS.vz - 1.0)*KS.fy/72. - (3.0*KS.vx - 3.0*KS.vy - 2.0*KS.vz + 1.0)*KS.fz/72.;
		const dreal Smzm = (2.0*KS.vx + 3.0*KS.vz - 1.0)*KS.fx/18. + (3.0*KS.vx + 2.0*KS.vz - 1.0)*KS.fz/18. - KS.fy*KS.vy/18.;
		const dreal Szzm = (4.0*KS.vz - 2.0)*KS.fz/9. - 2.0/9.0*(KS.fx*KS.vx + KS.vy*KS.fy);
		const dreal Spzm = (2.0*KS.vx - 3.0*KS.vz + 1.0)*KS.fx/18. + (-3.0*KS.vx + 2.0*KS.vz - 1.0)*KS.fz/18. - KS.fy*KS.vy/18.;
		const dreal Smpm = (2.0*KS.vx - 3.0*KS.vy + 3.0*KS.vz - 1.0)*KS.fx/72. + (-3.0*KS.vx + 2.0*KS.vy - 3.0*KS.vz + 1.0)*KS.fy/72. + (3.0*KS.vx - 3.0*KS.vy + 2.0*KS.vz - 1.0)*KS.fz/72.;
		const dreal Szpm = (2.0*KS.vy - 3.0*KS.vz + 1.0)*KS.fy/18. + (-3.0*KS.vy + 2.0*KS.vz - 1.0)*KS.fz/18. - KS.fx*KS.vx/18.;
		const dreal Sppm = (2.0*KS.vx + 3.0*KS.vy - 3.0*KS.vz + 1.0)*KS.fx/72. + (3.0*KS.vx + 2.0*KS.vy - 3.0*KS.vz + 1.0)*KS.fy/72. - (3.0*KS.vx + 3.0*KS.vy - 2.0*KS.vz + 1.0)*KS.fz/72.;
		const dreal Smmz = (2.0*KS.vx + 3.0*KS.vy - 1.0)*KS.fx/18. + (3.0*KS.vx + 2.0*KS.vy - 1.0)*KS.fy/18. - KS.vz*KS.fz/18.;
		const dreal Szmz = (4.0*KS.vy - 2.0)*KS.fy/9. - 2.0/9.0*(KS.vx*KS.fx + KS.vz*KS.fz);
		const dreal Spmz = (2.0*KS.vx - 3.0*KS.vy + 1.0)*KS.fx/18. + (-3.0*KS.vx + 2.0*KS.vy - 1.0)*KS.fy/18. - KS.fz*KS.vz/18.;
		const dreal Smzz = (4*KS.vx -2.0)*KS.fx/9.0 - 2.0/9.0*(KS.vy*KS.fy + KS.vz*KS.fz); 
		const dreal Szzz = -8.0/9.0*(KS.vx*KS.fx + KS.fy*KS.vy + KS.vz*KS.fz);
		const dreal Spzz = (4.0*KS.vx + 2.0)*KS.fx/9. - 2.0/9.0*(KS.vy*KS.fy + KS.vz*KS.fz);
		const dreal Smpz = (2.0*KS.vx - 3.0*KS.vy - 1.0)*KS.fx/18. + (-3.0*KS.vx + 2.0*KS.vy + 1.0)*KS.fy/18. - KS.fz*KS.vz/18.;
		const dreal Szpz = (4.0*KS.vy + 2.0)*KS.fy/9.0 - 2.0/9.0*(KS.fx*KS.vx + KS.fz*KS.vz);
		const dreal Sppz = (2.0*KS.vx + 3.0*KS.vy + 1.0)*KS.fx/18. + (3.0*KS.vx + 2.0*KS.vy + 1.0)*KS.fy/18. - KS.vz*KS.fz/18.;
		const dreal Smmp = (2.0*KS.vx + 3.0*KS.vy - 3.0*KS.vz - 1.0)*KS.fx/72. + (3.0*KS.vx + 2.0*KS.vy - 3.0*KS.vz - 1.0)*KS.fy/72. - (3.0*KS.vx + 3.0*KS.vy - 2.0*KS.vz - 1.0)*KS.fz/72.;
		const dreal Szmp = (2.*KS.vy - 3.*KS.vz - 1.)*KS.fy/18. + (-3.*KS.vy + 2.*KS.vz + 1.)*KS.fz/18. - KS.fx*KS.vx/18.; 
		const dreal Spmp = (2.0*KS.vx - 3.0*KS.vy + 3.0*KS.vz + 1.0)*KS.fx/72. + (-3.0*KS.vx + 2.0*KS.vy - 3.0*KS.vz - 1.0)*KS.fy/72. + (3.0*KS.vx - 3.0*KS.vy + 2.0*KS.vz + 1.0)*KS.fz/72.;
		const dreal Smzp = (2.0*KS.vx - 3.0*KS.vz - 1.0)*KS.fx/18. + (-3.0*KS.vx + 2.0*KS.vz + 1.0)*KS.fz/18. - KS.fy*KS.vy/18.; 
		const dreal Szzp = (4.0*KS.vz + 2.0)*KS.fz/9.0 - 2.0/9.0*(KS.fx*KS.vx + KS.vy*KS.fy);
		const dreal Spzp = (2.0*KS.vx + 3.0*KS.vz + 1.0)*KS.fx/18. + (3.0*KS.vx + 2.0*KS.vz + 1.0)*KS.fz/18. - KS.fy*KS.vy/18.;
		const dreal Smpp = (2.0*KS.vx - 3.0*KS.vy - 3.0*KS.vz - 1.0)*KS.fx/72. + (-3.0*KS.vx + 2.0*KS.vy + 3.0*KS.vz + 1.0)*KS.fy/72. - (3.0*KS.vx - 3.0*KS.vy - 2.0*KS.vz - 1.0)*KS.fz/72.;
		const dreal Szpp = (2.0*KS.vy + 3.0*KS.vz + 1.0)*KS.fy/18. + (3.0*KS.vy + 2.0*KS.vz + 1.0)*KS.fz/18. - KS.vx*KS.fx/18.;
		const dreal Sppp = (2.0*KS.vx + 3.0*KS.vy + 3.0*KS.vz + 1.0)*KS.fx/72. + (3.0*KS.vx + 2.0*KS.vy + 3.0*KS.vz + 1.0)*KS.fy/72. + (3.0*KS.vx + 3.0*KS.vy + 2.0*KS.vz + 1.0)*KS.fz/72.;

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

		KS.f[mmm] += (locfeq_mmm - KS.f[mmm])/tau + (no1 - n1o2/tau)*Smmm;
		KS.f[mmz] += (locfeq_mmz - KS.f[mmz])/tau + (no1 - n1o2/tau)*Smmz;
		KS.f[mmp] += (locfeq_mmp - KS.f[mmp])/tau + (no1 - n1o2/tau)*Smmp;
		KS.f[mzm] += (locfeq_mzm - KS.f[mzm])/tau + (no1 - n1o2/tau)*Smzm;
		KS.f[mzz] += (locfeq_mzz - KS.f[mzz])/tau + (no1 - n1o2/tau)*Smzz;
		KS.f[mzp] += (locfeq_mzp - KS.f[mzp])/tau + (no1 - n1o2/tau)*Smzp;
		KS.f[mpm] += (locfeq_mpm - KS.f[mpm])/tau + (no1 - n1o2/tau)*Smpm;
		KS.f[mpz] += (locfeq_mpz - KS.f[mpz])/tau + (no1 - n1o2/tau)*Smpz;
		KS.f[mpp] += (locfeq_mpp - KS.f[mpp])/tau + (no1 - n1o2/tau)*Smpp;
		KS.f[zmm] += (locfeq_zmm - KS.f[zmm])/tau + (no1 - n1o2/tau)*Szmm;
		KS.f[zmz] += (locfeq_zmz - KS.f[zmz])/tau + (no1 - n1o2/tau)*Szmz;
		KS.f[zmp] += (locfeq_zmp - KS.f[zmp])/tau + (no1 - n1o2/tau)*Szmp;
		KS.f[zzm] += (locfeq_zzm - KS.f[zzm])/tau + (no1 - n1o2/tau)*Szzm;
		KS.f[zzz] += (locfeq_zzz - KS.f[zzz])/tau + (no1 - n1o2/tau)*Szzz;
		KS.f[zzp] += (locfeq_zzp - KS.f[zzp])/tau + (no1 - n1o2/tau)*Szzp;
		KS.f[zpm] += (locfeq_zpm - KS.f[zpm])/tau + (no1 - n1o2/tau)*Szpm;
		KS.f[zpz] += (locfeq_zpz - KS.f[zpz])/tau + (no1 - n1o2/tau)*Szpz;
		KS.f[zpp] += (locfeq_zpp - KS.f[zpp])/tau + (no1 - n1o2/tau)*Szpp;
		KS.f[pmm] += (locfeq_pmm - KS.f[pmm])/tau + (no1 - n1o2/tau)*Spmm;
		KS.f[pmz] += (locfeq_pmz - KS.f[pmz])/tau + (no1 - n1o2/tau)*Spmz;
		KS.f[pmp] += (locfeq_pmp - KS.f[pmp])/tau + (no1 - n1o2/tau)*Spmp;
		KS.f[pzm] += (locfeq_pzm - KS.f[pzm])/tau + (no1 - n1o2/tau)*Spzm;
		KS.f[pzz] += (locfeq_pzz - KS.f[pzz])/tau + (no1 - n1o2/tau)*Spzz;
		KS.f[pzp] += (locfeq_pzp - KS.f[pzp])/tau + (no1 - n1o2/tau)*Spzp;
		KS.f[ppm] += (locfeq_ppm - KS.f[ppm])/tau + (no1 - n1o2/tau)*Sppm;
		KS.f[ppz] += (locfeq_ppz - KS.f[ppz])/tau + (no1 - n1o2/tau)*Sppz;
		KS.f[ppp] += (locfeq_ppp - KS.f[ppp])/tau + (no1 - n1o2/tau)*Sppp;
	}
};
