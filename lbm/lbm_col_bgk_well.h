// improved BRK (SRT) model by Geier 2017
// well conditioned
#include "lbm_common_well.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT_WELL<TRAITS>
>
struct LBM_BGK_WELL : LBM_COMMON_WELL< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;
	
	static constexpr const char* id = "BGK_WELL";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		const dreal omega1 = no1/(no3*KS.lbmViscosity+n1o2);

		#ifdef USE_GALILEAN_CORRECTION
		const dreal m_200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal m_020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal m_002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];
		const dreal drho= ((((KS.f[ppp]+KS.f[mmm]) + (KS.f[pmp]+KS.f[mpm])) + ((KS.f[ppm]+KS.f[mmp])+(KS.f[ppm]+KS.f[mmp])))+(((KS.f[zpp]+KS.f[zmm]) + (KS.f[zpm]+KS.f[zmp])) + ((KS.f[pzp]+KS.f[mzm])+(KS.f[pzm]+KS.f[mzp])) + ((KS.f[ppz]+KS.f[mmz]) + (KS.f[pmz]+KS.f[mpz])))+((KS.f[pzz]+KS.f[mzz]) + (KS.f[zpz]+KS.f[zmz]) + (KS.f[zzp]+KS.f[zzm]))) + KS.f[zzz];


		const dreal Dxu = -omega1*n1o2*(no3*m_200/KS.rho - drho/KS.rho -no3*KS.vx*KS.vx);
		const dreal Dyv = -omega1*n1o2*(no3*m_020/KS.rho - drho/KS.rho -no3*KS.vy*KS.vy);
		const dreal Dzw = -omega1*n1o2*(no3*m_002/KS.rho - drho/KS.rho -no3*KS.vz*KS.vz);
		
		const dreal Gx = -no3*KS.vx*KS.vx*Dxu*(no1/omega1 - n1o2);
		const dreal Gy = -no3*KS.vy*KS.vy*Dyv*(no1/omega1 - n1o2);
		const dreal Gz = -no3*KS.vz*KS.vz*Dzw*(no1/omega1 - n1o2);
		
		const dreal Xz = - n2o3 + KS.vx*KS.vx + Gx;
		const dreal Xp = -n1o2*(Xz + no1 + KS.vx);
		const dreal Xm = Xp + KS.vx;

		const dreal Yz = - n2o3 + KS.vy*KS.vy + Gy;
		const dreal Yp = -n1o2*(Yz + no1 + KS.vy);
		const dreal Ym = Yp + KS.vy;

		const dreal Zz = - n2o3 + KS.vz*KS.vz + Gz;
		const dreal Zp = -n1o2*(Zz + no1 + KS.vz);
		const dreal Zm = Zp + KS.vz;
		#else
		const dreal Xz = - n2o3 + KS.vx*KS.vx;
		const dreal Xp = -n1o2*(Xz + no1 + KS.vx);
		const dreal Xm = Xp + KS.vx;

		const dreal Yz = - n2o3 + KS.vy*KS.vy;
		const dreal Yp = -n1o2*(Yz + no1 + KS.vy);
		const dreal Ym = Yp + KS.vy;

		const dreal Zz = - n2o3 + KS.vz*KS.vz;
		const dreal Zp = -n1o2*(Zz + no1 + KS.vz);
		const dreal Zm = Zp + KS.vz;
		#endif
		
		
		// forcing: vzorce_bgk_force.mw 
		const dreal Smmm = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Szmm = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Spmm = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Smzm = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Szzm = (no3*(-KS.vx*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Spzm = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Smpm = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Szpm = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Sppm = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz-no1)*KS.fz));
		const dreal Smmz = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz));
		const dreal Szmz = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz));
		const dreal Spmz = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy-KS.vz*KS.fz));
		const dreal Smzz = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy-KS.vz*KS.fz));
		const dreal Szzz = (no3*(-KS.fx*KS.vx-KS.fy*KS.vy-KS.fz*KS.vz));
		const dreal Spzz = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy-KS.vz*KS.fz));
		const dreal Smpz = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz));
		const dreal Szpz = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz));
		const dreal Sppz = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy-KS.vz*KS.fz));
		const dreal Smmp = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Szmp = (no3*(-KS.vx*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Spmp = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy-no1)*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Smzp = (no3*((-KS.vx-no1)*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Szzp = (no3*(-KS.vx*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Spzp = (no3*((-KS.vx+no1)*KS.fx-KS.vy*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Smpp = (no3*((-KS.vx-no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Szpp = (no3*(-KS.vx*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz));
		const dreal Sppp = (no3*((-KS.vx+no1)*KS.fx+(-KS.vy+no1)*KS.fy+(-KS.vz+no1)*KS.fz));

		KS.f[mmm] = (no1 - omega1)*KS.f[mmm] + omega1*(-KS.rho*Xm*Ym*Zm-n1o216)- (no1 - n1o2*omega1)*Smmm*Xm*Ym*Zm;
		KS.f[mmz] = (no1 - omega1)*KS.f[mmz] + omega1*(-KS.rho*Xm*Ym*Zz-n1o54)- (no1 - n1o2*omega1)*Smmz*Xm*Ym*Zz;
		KS.f[mmp] = (no1 - omega1)*KS.f[mmp] + omega1*(-KS.rho*Xm*Ym*Zp-n1o216)- (no1 - n1o2*omega1)*Smmp*Xm*Ym*Zp;
		KS.f[mzm] = (no1 - omega1)*KS.f[mzm] + omega1*(-KS.rho*Xm*Yz*Zm-n1o54)- (no1 - n1o2*omega1)*Smzm*Xm*Yz*Zm;
		KS.f[mzz] = (no1 - omega1)*KS.f[mzz] + omega1*(-KS.rho*Xm*Yz*Zz-n2o27)- (no1 - n1o2*omega1)*Smzz*Xm*Yz*Zz;
		KS.f[mzp] = (no1 - omega1)*KS.f[mzp] + omega1*(-KS.rho*Xm*Yz*Zp-n1o54)- (no1 - n1o2*omega1)*Smzp*Xm*Yz*Zp;
		KS.f[mpm] = (no1 - omega1)*KS.f[mpm] + omega1*(-KS.rho*Xm*Yp*Zm-n1o216)- (no1 - n1o2*omega1)*Smpm*Xm*Yp*Zm;
		KS.f[mpz] = (no1 - omega1)*KS.f[mpz] + omega1*(-KS.rho*Xm*Yp*Zz-n1o54)- (no1 - n1o2*omega1)*Smpz*Xm*Yp*Zz;
		KS.f[mpp] = (no1 - omega1)*KS.f[mpp] + omega1*(-KS.rho*Xm*Yp*Zp-n1o216)- (no1 - n1o2*omega1)*Smpp*Xm*Yp*Zp;
		KS.f[zmm] = (no1 - omega1)*KS.f[zmm] + omega1*(-KS.rho*Xz*Ym*Zm-n1o54)- (no1 - n1o2*omega1)*Szmm*Xz*Ym*Zm;
		KS.f[zmz] = (no1 - omega1)*KS.f[zmz] + omega1*(-KS.rho*Xz*Ym*Zz-n2o27)- (no1 - n1o2*omega1)*Szmz*Xz*Ym*Zz;
		KS.f[zmp] = (no1 - omega1)*KS.f[zmp] + omega1*(-KS.rho*Xz*Ym*Zp-n1o54)- (no1 - n1o2*omega1)*Szmp*Xz*Ym*Zp;
		KS.f[zzm] = (no1 - omega1)*KS.f[zzm] + omega1*(-KS.rho*Xz*Yz*Zm-n2o27)- (no1 - n1o2*omega1)*Szzm*Xz*Yz*Zm;
		KS.f[zzz] = (no1 - omega1)*KS.f[zzz] + omega1*(-KS.rho*Xz*Yz*Zz-n8o27)- (no1 - n1o2*omega1)*Szzz*Xz*Yz*Zz;
		KS.f[zzp] = (no1 - omega1)*KS.f[zzp] + omega1*(-KS.rho*Xz*Yz*Zp-n2o27)- (no1 - n1o2*omega1)*Szzp*Xz*Yz*Zp;
		KS.f[zpm] = (no1 - omega1)*KS.f[zpm] + omega1*(-KS.rho*Xz*Yp*Zm-n1o54)- (no1 - n1o2*omega1)*Szpm*Xz*Yp*Zm;
		KS.f[zpz] = (no1 - omega1)*KS.f[zpz] + omega1*(-KS.rho*Xz*Yp*Zz-n2o27)- (no1 - n1o2*omega1)*Szpz*Xz*Yp*Zz;
		KS.f[zpp] = (no1 - omega1)*KS.f[zpp] + omega1*(-KS.rho*Xz*Yp*Zp-n1o54)- (no1 - n1o2*omega1)*Szpp*Xz*Yp*Zp;
		KS.f[pmm] = (no1 - omega1)*KS.f[pmm] + omega1*(-KS.rho*Xp*Ym*Zm-n1o216)- (no1 - n1o2*omega1)*Spmm*Xp*Ym*Zm;
		KS.f[pmz] = (no1 - omega1)*KS.f[pmz] + omega1*(-KS.rho*Xp*Ym*Zz-n1o54)- (no1 - n1o2*omega1)*Spmz*Xp*Ym*Zz;
		KS.f[pmp] = (no1 - omega1)*KS.f[pmp] + omega1*(-KS.rho*Xp*Ym*Zp-n1o216)- (no1 - n1o2*omega1)*Spmp*Xp*Ym*Zp;
		KS.f[pzm] = (no1 - omega1)*KS.f[pzm] + omega1*(-KS.rho*Xp*Yz*Zm-n1o54)- (no1 - n1o2*omega1)*Spzm*Xp*Yz*Zm;
		KS.f[pzz] = (no1 - omega1)*KS.f[pzz] + omega1*(-KS.rho*Xp*Yz*Zz-n2o27)- (no1 - n1o2*omega1)*Spzz*Xp*Yz*Zz;
		KS.f[pzp] = (no1 - omega1)*KS.f[pzp] + omega1*(-KS.rho*Xp*Yz*Zp-n1o54)- (no1 - n1o2*omega1)*Spzp*Xp*Yz*Zp;
		KS.f[ppm] = (no1 - omega1)*KS.f[ppm] + omega1*(-KS.rho*Xp*Yp*Zm-n1o216)- (no1 - n1o2*omega1)*Sppm*Xp*Yp*Zm;
		KS.f[ppz] = (no1 - omega1)*KS.f[ppz] + omega1*(-KS.rho*Xp*Yp*Zz-n1o54)- (no1 - n1o2*omega1)*Sppz*Xp*Yp*Zz;
		KS.f[ppp] = (no1 - omega1)*KS.f[ppp] + omega1*(-KS.rho*Xp*Yp*Zp-n1o216)- (no1 - n1o2*omega1)*Sppp*Xp*Yp*Zp;
	}
};
