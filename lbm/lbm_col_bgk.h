// improved BRK (SRT) model by Geier 2017
// for standard DF (no well-conditioned)
#include "lbm_common.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT< TRAITS >
>
struct LBM_BGK : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "BGK";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		const dreal omega1 = no1/(no3*KS.lbmViscosity+n1o2);
		#ifdef USE_GALILEAN_CORRECTION
		const dreal m_200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal m_020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal m_002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];


		const dreal Dxu = -omega1*n1o2*(no3*m_200/KS.rho - no1 -no3*KS.vx*KS.vx);
		const dreal Dyv = -omega1*n1o2*(no3*m_020/KS.rho - no1 -no3*KS.vy*KS.vy);
		const dreal Dzw = -omega1*n1o2*(no3*m_002/KS.rho - no1 -no3*KS.vz*KS.vz);

		const dreal Gx = -no3*KS.vx*KS.vx*Dxu*(no1/omega1 - n1o2);
		const dreal Gy = -no3*KS.vy*KS.vy*Dyv*(no1/omega1 - n1o2);
		const dreal Gz = -no3*KS.vz*KS.vz*Dzw*(no1/omega1 - n1o2);
		
		const dreal Xz = n1o3 - no1 + KS.vx*KS.vx + Gx;
		const dreal Xp = -n1o2*(Xz + no1 + KS.vx);
		const dreal Xm = Xp + KS.vx;

		const dreal Yz = n1o3 - no1 + KS.vy*KS.vy + Gy;
		const dreal Yp = -n1o2*(Yz + no1 + KS.vy);
		const dreal Ym = Yp + KS.vy;

		const dreal Zz = n1o3 - no1 + KS.vz*KS.vz + Gz;
		const dreal Zp = -n1o2*(Zz + no1 + KS.vz);
		const dreal Zm = Zp + KS.vz;
		#else
		const dreal Xz = n1o3 - no1 + KS.vx*KS.vx;
		const dreal Xp = -n1o2*(Xz + no1 + KS.vx);
		const dreal Xm = Xp + KS.vx;

		const dreal Yz = n1o3 - no1 + KS.vy*KS.vy;
		const dreal Yp = -n1o2*(Yz + no1 + KS.vy);
		const dreal Ym = Yp + KS.vy;

		const dreal Zz = n1o3 - no1 + KS.vz*KS.vz;
		const dreal Zp = -n1o2*(Zz + no1 + KS.vz);
		const dreal Zm = Zp + KS.vz;
		#endif

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

		const dreal feq_mmm = -KS.rho*Xm*Ym*Zm;
		const dreal feq_mmz = -KS.rho*Xm*Ym*Zz;
		const dreal feq_mmp = -KS.rho*Xm*Ym*Zp;
		const dreal feq_mzm = -KS.rho*Xm*Yz*Zm;
		const dreal feq_mzz = -KS.rho*Xm*Yz*Zz;
		const dreal feq_mzp = -KS.rho*Xm*Yz*Zp;
		const dreal feq_mpm = -KS.rho*Xm*Yp*Zm;
		const dreal feq_mpz = -KS.rho*Xm*Yp*Zz;
		const dreal feq_mpp = -KS.rho*Xm*Yp*Zp;
		const dreal feq_zmm = -KS.rho*Xz*Ym*Zm;
		const dreal feq_zmz = -KS.rho*Xz*Ym*Zz;
		const dreal feq_zmp = -KS.rho*Xz*Ym*Zp;
		const dreal feq_zzm = -KS.rho*Xz*Yz*Zm;
		const dreal feq_zzz = -KS.rho*Xz*Yz*Zz;
		const dreal feq_zzp = -KS.rho*Xz*Yz*Zp;
		const dreal feq_zpm = -KS.rho*Xz*Yp*Zm;
		const dreal feq_zpz = -KS.rho*Xz*Yp*Zz;
		const dreal feq_zpp = -KS.rho*Xz*Yp*Zp;
		const dreal feq_pmm = -KS.rho*Xp*Ym*Zm;
		const dreal feq_pmz = -KS.rho*Xp*Ym*Zz;
		const dreal feq_pmp = -KS.rho*Xp*Ym*Zp;
		const dreal feq_pzm = -KS.rho*Xp*Yz*Zm;
		const dreal feq_pzz = -KS.rho*Xp*Yz*Zz;
		const dreal feq_pzp = -KS.rho*Xp*Yz*Zp;
		const dreal feq_ppm = -KS.rho*Xp*Yp*Zm;
		const dreal feq_ppz = -KS.rho*Xp*Yp*Zz;
		const dreal feq_ppp = -KS.rho*Xp*Yp*Zp;

		KS.f[mmm] += (feq_mmm - KS.f[mmm])*omega1 + (no1 - n1o2*omega1)*Smmm*feq_mmm;
		KS.f[mmz] += (feq_mmz - KS.f[mmz])*omega1 + (no1 - n1o2*omega1)*Smmz*feq_mmz;
		KS.f[mmp] += (feq_mmp - KS.f[mmp])*omega1 + (no1 - n1o2*omega1)*Smmp*feq_mmp;
		KS.f[mzm] += (feq_mzm - KS.f[mzm])*omega1 + (no1 - n1o2*omega1)*Smzm*feq_mzm;
		KS.f[mzz] += (feq_mzz - KS.f[mzz])*omega1 + (no1 - n1o2*omega1)*Smzz*feq_mzz;
		KS.f[mzp] += (feq_mzp - KS.f[mzp])*omega1 + (no1 - n1o2*omega1)*Smzp*feq_mzp;
		KS.f[mpm] += (feq_mpm - KS.f[mpm])*omega1 + (no1 - n1o2*omega1)*Smpm*feq_mpm;
		KS.f[mpz] += (feq_mpz - KS.f[mpz])*omega1 + (no1 - n1o2*omega1)*Smpz*feq_mpz;
		KS.f[mpp] += (feq_mpp - KS.f[mpp])*omega1 + (no1 - n1o2*omega1)*Smpp*feq_mpp;
		KS.f[zmm] += (feq_zmm - KS.f[zmm])*omega1 + (no1 - n1o2*omega1)*Szmm*feq_zmm;
		KS.f[zmz] += (feq_zmz - KS.f[zmz])*omega1 + (no1 - n1o2*omega1)*Szmz*feq_zmz;
		KS.f[zmp] += (feq_zmp - KS.f[zmp])*omega1 + (no1 - n1o2*omega1)*Szmp*feq_zmp;
		KS.f[zzm] += (feq_zzm - KS.f[zzm])*omega1 + (no1 - n1o2*omega1)*Szzm*feq_zzm;
		KS.f[zzz] += (feq_zzz - KS.f[zzz])*omega1 + (no1 - n1o2*omega1)*Szzz*feq_zzz;
		KS.f[zzp] += (feq_zzp - KS.f[zzp])*omega1 + (no1 - n1o2*omega1)*Szzp*feq_zzp;
		KS.f[zpm] += (feq_zpm - KS.f[zpm])*omega1 + (no1 - n1o2*omega1)*Szpm*feq_zpm;
		KS.f[zpz] += (feq_zpz - KS.f[zpz])*omega1 + (no1 - n1o2*omega1)*Szpz*feq_zpz;
		KS.f[zpp] += (feq_zpp - KS.f[zpp])*omega1 + (no1 - n1o2*omega1)*Szpp*feq_zpp;
		KS.f[pmm] += (feq_pmm - KS.f[pmm])*omega1 + (no1 - n1o2*omega1)*Spmm*feq_pmm;
		KS.f[pmz] += (feq_pmz - KS.f[pmz])*omega1 + (no1 - n1o2*omega1)*Spmz*feq_pmz;
		KS.f[pmp] += (feq_pmp - KS.f[pmp])*omega1 + (no1 - n1o2*omega1)*Spmp*feq_pmp;
		KS.f[pzm] += (feq_pzm - KS.f[pzm])*omega1 + (no1 - n1o2*omega1)*Spzm*feq_pzm;
		KS.f[pzz] += (feq_pzz - KS.f[pzz])*omega1 + (no1 - n1o2*omega1)*Spzz*feq_pzz;
		KS.f[pzp] += (feq_pzp - KS.f[pzp])*omega1 + (no1 - n1o2*omega1)*Spzp*feq_pzp;
		KS.f[ppm] += (feq_ppm - KS.f[ppm])*omega1 + (no1 - n1o2*omega1)*Sppm*feq_ppm;
		KS.f[ppz] += (feq_ppz - KS.f[ppz])*omega1 + (no1 - n1o2*omega1)*Sppz*feq_ppz;
		KS.f[ppp] += (feq_ppp - KS.f[ppp])*omega1 + (no1 - n1o2*omega1)*Sppp*feq_ppp;

/*		
		KS.f[mmm] = KS.f[mmm]*(no1-omega1) + omega1*(-KS.rho*Xm*Ym*Zm);
		KS.f[mmz] = KS.f[mmz]*(no1-omega1) + omega1*(-KS.rho*Xm*Ym*Zz);
		KS.f[mmp] = KS.f[mmp]*(no1-omega1) + omega1*(-KS.rho*Xm*Ym*Zp);
		KS.f[mzm] = KS.f[mzm]*(no1-omega1) + omega1*(-KS.rho*Xm*Yz*Zm);
		KS.f[mzz] = KS.f[mzz]*(no1-omega1) + omega1*(-KS.rho*Xm*Yz*Zz);
		KS.f[mzp] = KS.f[mzp]*(no1-omega1) + omega1*(-KS.rho*Xm*Yz*Zp);
		KS.f[mpm] = KS.f[mpm]*(no1-omega1) + omega1*(-KS.rho*Xm*Yp*Zm);
		KS.f[mpz] = KS.f[mpz]*(no1-omega1) + omega1*(-KS.rho*Xm*Yp*Zz);
		KS.f[mpp] = KS.f[mpp]*(no1-omega1) + omega1*(-KS.rho*Xm*Yp*Zp);
		KS.f[zmm] = KS.f[zmm]*(no1-omega1) + omega1*(-KS.rho*Xz*Ym*Zm);
		KS.f[zmz] = KS.f[zmz]*(no1-omega1) + omega1*(-KS.rho*Xz*Ym*Zz);
		KS.f[zmp] = KS.f[zmp]*(no1-omega1) + omega1*(-KS.rho*Xz*Ym*Zp);
		KS.f[zzm] = KS.f[zzm]*(no1-omega1) + omega1*(-KS.rho*Xz*Yz*Zm);
		KS.f[zzz] = KS.f[zzz]*(no1-omega1) + omega1*(-KS.rho*Xz*Yz*Zz);
		KS.f[zzp] = KS.f[zzp]*(no1-omega1) + omega1*(-KS.rho*Xz*Yz*Zp);
		KS.f[zpm] = KS.f[zpm]*(no1-omega1) + omega1*(-KS.rho*Xz*Yp*Zm);
		KS.f[zpz] = KS.f[zpz]*(no1-omega1) + omega1*(-KS.rho*Xz*Yp*Zz);
		KS.f[zpp] = KS.f[zpp]*(no1-omega1) + omega1*(-KS.rho*Xz*Yp*Zp);
		KS.f[pmm] = KS.f[pmm]*(no1-omega1) + omega1*(-KS.rho*Xp*Ym*Zm);
		KS.f[pmz] = KS.f[pmz]*(no1-omega1) + omega1*(-KS.rho*Xp*Ym*Zz);
		KS.f[pmp] = KS.f[pmp]*(no1-omega1) + omega1*(-KS.rho*Xp*Ym*Zp);
		KS.f[pzm] = KS.f[pzm]*(no1-omega1) + omega1*(-KS.rho*Xp*Yz*Zm);
		KS.f[pzz] = KS.f[pzz]*(no1-omega1) + omega1*(-KS.rho*Xp*Yz*Zz);
		KS.f[pzp] = KS.f[pzp]*(no1-omega1) + omega1*(-KS.rho*Xp*Yz*Zp);
		KS.f[ppm] = KS.f[ppm]*(no1-omega1) + omega1*(-KS.rho*Xp*Yp*Zm);
		KS.f[ppz] = KS.f[ppz]*(no1-omega1) + omega1*(-KS.rho*Xp*Yp*Zz);
		KS.f[ppp] = KS.f[ppp]*(no1-omega1) + omega1*(-KS.rho*Xp*Yp*Zp);
*/
	}
};
