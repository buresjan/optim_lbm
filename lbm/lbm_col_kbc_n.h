// KBC models by Karlin Bosch Chikatamarla 2015
/* arXiv:1507.02518v1 [physics.flu-dyn] 

 * for standard DF (no well-conditioned) D3Q27 with c_s^2 = 1/3
 * eight models can be derived with different shear part for the DF decomposition f_i = k_i + s_i + h_i
 * N1,N2,N3,N4 with natural/raw moments and C1,C2,C3,C4 with central moment representation of DF
 * bulk visc. is equal to kin. for N2,N4,C2,C4 (i.e. when shear part contains trace of pressure tensor)
 * or equal to cs^2(1/(gamma*beta)-.5) for N1,N3,C1,C3
 * shear part s_i consists of following tensors:
 * N1/C1 D/D~
 * N2/C2 D+T/D~+T~
 * N3/C3 D+Q/D~+Q~
 * N4/C4 D+T+Q/D~+T~+Q~
 * see the above paper (and few others by K+B+C) for more details on D & T (and not only!), tilde ~ means central moments
 */
// N1-N4 models with the raw/natural moment representation of DFs

#include "lbm_common.h"

// speed of sound
#define _csqr n1o3

//raw moments - models N1, N2, N3, N4
// special moments with physical meaning & their equilibria
#define rT    (M200 + M020 + M002)
#define rNxz  (M200 - M002)
#define rNyz  (M020 - M002)
#define rPxy  M110
#define rPxz  M101
#define rPyz  M011
#define rQxxy M210
#define rQxxz M201
#define rQxyy M120
#define rQyyz M021
#define rQxzz M102
#define rQyzz M012
#define rQxyz M111
#define erT    (KS.rho*(no3*_csqr + KS.vx*KS.vx + KS.vy*KS.vy + KS.vz*KS.vz))
#define erNxz  (KS.rho*(KS.vx*KS.vx - KS.vz*KS.vz))
#define erNyz  (KS.rho*(KS.vy*KS.vy - KS.vz*KS.vz))
#define erPxy  (KS.rho*KS.vx*KS.vy)
#define erPxz  (KS.rho*KS.vx*KS.vz)
#define erPyz  (KS.rho*KS.vy*KS.vz)
#define erQxxy (KS.rho*KS.vy*(_csqr + KS.vx*KS.vx))
#define erQxxz (KS.rho*KS.vz*(_csqr + KS.vx*KS.vx))
#define erQxyy (KS.rho*KS.vx*(_csqr + KS.vy*KS.vy))
#define erQyyz (KS.rho*KS.vz*(_csqr + KS.vy*KS.vy))
#define erQxzz (KS.rho*KS.vx*(_csqr + KS.vz*KS.vz))
#define erQyzz (KS.rho*KS.vy*(_csqr + KS.vz*KS.vz))
#define erQxyz (KS.rho*KS.vx*KS.vy*KS.vz)
// tensors included in the shear part
//  deviatoric stress tensor D
#define Dzzz no0
#define Dpzz ((no2*rNxz - rNyz)*n1o6)
#define Dmzz ((no2*rNxz - rNyz)*n1o6)
#define Dzpz ((-rNxz + no2*rNyz)*n1o6)
#define Dzmz ((-rNxz + no2*rNyz)*n1o6)
#define Dzzp ((-rNxz - rNyz)*n1o6)
#define Dzzm ((-rNxz - rNyz)*n1o6)
#define Dmmz (rPxy*n1o4)
#define Dmpz (-rPxy*n1o4)
#define Dpmz (-rPxy*n1o4)
#define Dppz (rPxy*n1o4)
#define Dmzm (rPxz*n1o4)
#define Dmzp (-rPxz*n1o4)
#define Dpzm (-rPxz*n1o4)
#define Dpzp (rPxz*n1o4)
#define Dzmm (rPyz*n1o4)
#define Dzmp (-rPyz*n1o4)
#define Dzpm (-rPyz*n1o4)
#define Dzpp (rPyz*n1o4)
#define Dmmm no0
#define Dmmp no0
#define Dmpm no0
#define Dmpp no0
#define Dpmm no0
#define Dpmp no0
#define Dppm no0
#define Dppp no0

#define eDzzz no0
#define eDpzz ((no2*erNxz - erNyz)*n1o6)
#define eDmzz ((no2*erNxz - erNyz)*n1o6)
#define eDzpz ((-erNxz + no2*erNyz)*n1o6)
#define eDzmz ((-erNxz + no2*erNyz)*n1o6)
#define eDzzp ((-erNxz - erNyz)*n1o6)
#define eDzzm ((-erNxz - erNyz)*n1o6)
#define eDmmz (erPxy*n1o4)
#define eDmpz (-erPxy*n1o4)
#define eDpmz (-erPxy*n1o4)
#define eDppz (erPxy*n1o4)
#define eDmzm (erPxz*n1o4)
#define eDmzp (-erPxz*n1o4)
#define eDpzm (-erPxz*n1o4)
#define eDpzp (erPxz*n1o4)
#define eDzmm (erPyz*n1o4)
#define eDzmp (-erPyz*n1o4)
#define eDzpm (-erPyz*n1o4)
#define eDzpp (erPyz*n1o4)
#define eDmmm no0
#define eDmmp no0
#define eDmpm no0
#define eDmpp no0
#define eDpmm no0
#define eDpmp no0
#define eDppm no0
#define eDppp no0
//  trace of stress tensor T
#define Tzzz -rT
#define Tpzz (rT*n1o6)
#define Tmzz (rT*n1o6)
#define Tzpz (rT*n1o6)
#define Tzmz (rT*n1o6)
#define Tzzp (rT*n1o6)
#define Tzzm (rT*n1o6)
#define Tmmz no0
#define Tmpz no0
#define Tpmz no0
#define Tppz no0
#define Tmzm no0
#define Tmzp no0
#define Tpzm no0
#define Tpzp no0
#define Tzmm no0
#define Tzmp no0
#define Tzpm no0
#define Tzpp no0
#define Tmmm no0
#define Tmmp no0
#define Tmpm no0
#define Tmpp no0
#define Tpmm no0
#define Tpmp no0
#define Tppm no0
#define Tppp no0

#define eTzzz -erT
#define eTpzz (erT*n1o6)
#define eTmzz (erT*n1o6)
#define eTzpz (erT*n1o6)
#define eTzmz (erT*n1o6)
#define eTzzp (erT*n1o6)
#define eTzzm (erT*n1o6)
#define eTmmz no0
#define eTmpz no0
#define eTpmz no0
#define eTppz no0
#define eTmzm no0
#define eTmzp no0
#define eTpzm no0
#define eTpzp no0
#define eTzmm no0
#define eTzmp no0
#define eTzpm no0
#define eTzpp no0
#define eTmmm no0
#define eTmmp no0
#define eTmpm no0
#define eTmpp no0
#define eTpmm no0
#define eTpmp no0
#define eTppm no0
#define eTppp no0
//  heat flux tensor Q
#define Qzzz no0
#define Qpzz (-(rQxyy + rQxzz)*n1o2)
#define Qmzz (+(rQxyy + rQxzz)*n1o2)
#define Qzpz (-(rQxxy + rQyzz)*n1o2)
#define Qzmz (+(rQxxy + rQyzz)*n1o2)
#define Qzzp (-(rQxxz + rQyyz)*n1o2)
#define Qzzm (+(rQxxz + rQyyz)*n1o2)
#define Qmmz ((-rQxyy - rQxxy)*n1o4)
#define Qmpz ((-rQxyy + rQxxy)*n1o4)
#define Qpmz ((+rQxyy - rQxxy)*n1o4)
#define Qppz ((+rQxyy + rQxxy)*n1o4)
#define Qmzm ((-rQxzz - rQxxz)*n1o4)
#define Qmzp ((-rQxzz + rQxxz)*n1o4)
#define Qpzm ((+rQxzz - rQxxz)*n1o4)
#define Qpzp ((+rQxzz + rQxxz)*n1o4)
#define Qzmm ((-rQyzz - rQyyz)*n1o4)
#define Qzmp ((-rQyzz + rQyyz)*n1o4)
#define Qzpm ((+rQyzz - rQyyz)*n1o4)
#define Qzpp ((+rQyzz + rQyyz)*n1o4)
#define Qmmm (-rQxyz*n1o8)
#define Qmmp (+rQxyz*n1o8)
#define Qmpm (+rQxyz*n1o8)
#define Qmpp (-rQxyz*n1o8)
#define Qpmm (+rQxyz*n1o8)
#define Qpmp (-rQxyz*n1o8)
#define Qppm (-rQxyz*n1o8)
#define Qppp (+rQxyz*n1o8)

#define eQzzz no0
#define eQpzz (-(erQxyy + erQxzz)*n1o2)
#define eQmzz (+(erQxyy + erQxzz)*n1o2)
#define eQzpz (-(erQxxy + erQyzz)*n1o2)
#define eQzmz (+(erQxxy + erQyzz)*n1o2)
#define eQzzp (-(erQxxz + erQyyz)*n1o2)
#define eQzzm (+(erQxxz + erQyyz)*n1o2)
#define eQmmz ((-erQxyy - erQxxy)*n1o4)
#define eQmpz ((-erQxyy + erQxxy)*n1o4)
#define eQpmz ((+erQxyy - erQxxy)*n1o4)
#define eQppz ((+erQxyy + erQxxy)*n1o4)
#define eQmzm ((-erQxzz - erQxxz)*n1o4)
#define eQmzp ((-erQxzz + erQxxz)*n1o4)
#define eQpzm ((+erQxzz - erQxxz)*n1o4)
#define eQpzp ((+erQxzz + erQxxz)*n1o4)
#define eQzmm ((-erQyzz - erQyyz)*n1o4)
#define eQzmp ((-erQyzz + erQyyz)*n1o4)
#define eQzpm ((+erQyzz - erQyyz)*n1o4)
#define eQzpp ((+erQyzz + erQyyz)*n1o4)
#define eQmmm (-erQxyz*n1o8)
#define eQmmp (+erQxyz*n1o8)
#define eQmpm (+erQxyz*n1o8)
#define eQmpp (-erQxyz*n1o8)
#define eQpmm (+erQxyz*n1o8)
#define eQpmp (-erQxyz*n1o8)
#define eQppm (-erQxyz*n1o8)
#define eQppp (+erQxyz*n1o8)

//relaxation parameter for shear part of DF - beta
// nu_lb = cs^2/2*(beta-1)

#define _beta (no1/(no2*KS.lbmViscosity/_csqr+no1))

//entropic stabilizer
// gamma = 1/beta - (2 - 1/beta)<Ds|Dh>/<Dh|Dh>
// <x|y> = sum(i) x_i*y_i/feq_i
// Regularized LB is achieved by setting gamma = 1/beta
//#define _gamma (no1/beta)
// SRT-BGK is achieved by setting gamma = 2
//#define _gamma no2


#define _gamma (no1/beta - (no2 - no1/beta)*\
	(Dsmmm*Dhmmm*ifeq_mmm +\
	 Dsmmz*Dhmmz*ifeq_mmz +\
	 Dsmmp*Dhmmp*ifeq_mmp +\
	 Dsmzm*Dhmzm*ifeq_mzm +\
	 Dsmzz*Dhmzz*ifeq_mzz +\
	 Dsmzp*Dhmzp*ifeq_mzp +\
	 Dsmpm*Dhmpm*ifeq_mpm +\
	 Dsmpz*Dhmpz*ifeq_mpz +\
	 Dsmpp*Dhmpp*ifeq_mpp +\
	 Dszmm*Dhzmm*ifeq_zmm +\
	 Dszmz*Dhzmz*ifeq_zmz +\
	 Dszmp*Dhzmp*ifeq_zmp +\
	 Dszzm*Dhzzm*ifeq_zzm +\
	 Dszzz*Dhzzz*ifeq_zzz +\
	 Dszzp*Dhzzp*ifeq_zzp +\
	 Dszpm*Dhzpm*ifeq_zpm +\
	 Dszpz*Dhzpz*ifeq_zpz +\
	 Dszpp*Dhzpp*ifeq_zpp +\
	 Dspmm*Dhpmm*ifeq_pmm +\
	 Dspmz*Dhpmz*ifeq_pmz +\
	 Dspmp*Dhpmp*ifeq_pmp +\
	 Dspzm*Dhpzm*ifeq_pzm +\
	 Dspzz*Dhpzz*ifeq_pzz +\
	 Dspzp*Dhpzp*ifeq_pzp +\
	 Dsppm*Dhppm*ifeq_ppm +\
	 Dsppz*Dhppz*ifeq_ppz +\
	 Dsppp*Dhppp*ifeq_ppp)/\
	(Dhmmm*Dhmmm*ifeq_mmm +\
	 Dhmmz*Dhmmz*ifeq_mmz +\
	 Dhmmp*Dhmmp*ifeq_mmp +\
	 Dhmzm*Dhmzm*ifeq_mzm +\
	 Dhmzz*Dhmzz*ifeq_mzz +\
	 Dhmzp*Dhmzp*ifeq_mzp +\
	 Dhmpm*Dhmpm*ifeq_mpm +\
	 Dhmpz*Dhmpz*ifeq_mpz +\
	 Dhmpp*Dhmpp*ifeq_mpp +\
	 Dhzmm*Dhzmm*ifeq_zmm +\
	 Dhzmz*Dhzmz*ifeq_zmz +\
	 Dhzmp*Dhzmp*ifeq_zmp +\
	 Dhzzm*Dhzzm*ifeq_zzm +\
	 Dhzzz*Dhzzz*ifeq_zzz +\
	 Dhzzp*Dhzzp*ifeq_zzp +\
	 Dhzpm*Dhzpm*ifeq_zpm +\
	 Dhzpz*Dhzpz*ifeq_zpz +\
	 Dhzpp*Dhzpp*ifeq_zpp +\
	 Dhpmm*Dhpmm*ifeq_pmm +\
	 Dhpmz*Dhpmz*ifeq_pmz +\
	 Dhpmp*Dhpmp*ifeq_pmp +\
	 Dhpzm*Dhpzm*ifeq_pzm +\
	 Dhpzz*Dhpzz*ifeq_pzz +\
	 Dhpzp*Dhpzp*ifeq_pzp +\
	 Dhppm*Dhppm*ifeq_ppm +\
	 Dhppz*Dhppz*ifeq_ppz +\
	 Dhppp*Dhppp*ifeq_ppp))


template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_KBC_N1 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_N1";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
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
		const dreal ffeq_mmm = -KS.rho*Xm*Ym*Zm;
		const dreal ffeq_mmz = -KS.rho*Xm*Ym*Zz;
		const dreal ffeq_mmp = -KS.rho*Xm*Ym*Zp;
		const dreal ffeq_mzm = -KS.rho*Xm*Yz*Zm;
		const dreal ffeq_mzz = -KS.rho*Xm*Yz*Zz;
		const dreal ffeq_mzp = -KS.rho*Xm*Yz*Zp;
		const dreal ffeq_mpm = -KS.rho*Xm*Yp*Zm;
		const dreal ffeq_mpz = -KS.rho*Xm*Yp*Zz;
		const dreal ffeq_mpp = -KS.rho*Xm*Yp*Zp;
		const dreal ffeq_zmm = -KS.rho*Xz*Ym*Zm;
		const dreal ffeq_zmz = -KS.rho*Xz*Ym*Zz;
		const dreal ffeq_zmp = -KS.rho*Xz*Ym*Zp;
		const dreal ffeq_zzm = -KS.rho*Xz*Yz*Zm;
		const dreal ffeq_zzz = -KS.rho*Xz*Yz*Zz;
		const dreal ffeq_zzp = -KS.rho*Xz*Yz*Zp;
		const dreal ffeq_zpm = -KS.rho*Xz*Yp*Zm;
		const dreal ffeq_zpz = -KS.rho*Xz*Yp*Zz;
		const dreal ffeq_zpp = -KS.rho*Xz*Yp*Zp;
		const dreal ffeq_pmm = -KS.rho*Xp*Ym*Zm;
		const dreal ffeq_pmz = -KS.rho*Xp*Ym*Zz;
		const dreal ffeq_pmp = -KS.rho*Xp*Ym*Zp;
		const dreal ffeq_pzm = -KS.rho*Xp*Yz*Zm;
		const dreal ffeq_pzz = -KS.rho*Xp*Yz*Zz;
		const dreal ffeq_pzp = -KS.rho*Xp*Yz*Zp;
		const dreal ffeq_ppm = -KS.rho*Xp*Yp*Zm;
		const dreal ffeq_ppz = -KS.rho*Xp*Yp*Zz;
		const dreal ffeq_ppp = -KS.rho*Xp*Yp*Zp;

		const dreal	ifeq_mmm = no1/ffeq_mmm;
		const dreal ifeq_mmz = no1/ffeq_mmz;
		const dreal ifeq_mmp = no1/ffeq_mmp;
		const dreal ifeq_mzm = no1/ffeq_mzm;
		const dreal ifeq_mzz = no1/ffeq_mzz;
		const dreal ifeq_mzp = no1/ffeq_mzp;
		const dreal ifeq_mpm = no1/ffeq_mpm;
		const dreal ifeq_mpz = no1/ffeq_mpz;
		const dreal ifeq_mpp = no1/ffeq_mpp;
		const dreal ifeq_zmm = no1/ffeq_zmm;
		const dreal ifeq_zmz = no1/ffeq_zmz;
		const dreal ifeq_zmp = no1/ffeq_zmp;
		const dreal ifeq_zzm = no1/ffeq_zzm;
		const dreal ifeq_zzz = no1/ffeq_zzz;
		const dreal ifeq_zzp = no1/ffeq_zzp;
		const dreal ifeq_zpm = no1/ffeq_zpm;
		const dreal ifeq_zpz = no1/ffeq_zpz;
		const dreal ifeq_zpp = no1/ffeq_zpp;
		const dreal ifeq_pmm = no1/ffeq_pmm;
		const dreal ifeq_pmz = no1/ffeq_pmz;
		const dreal ifeq_pmp = no1/ffeq_pmp;
		const dreal ifeq_pzm = no1/ffeq_pzm;
		const dreal ifeq_pzz = no1/ffeq_pzz;
		const dreal ifeq_pzp = no1/ffeq_pzp;
		const dreal ifeq_ppm = no1/ffeq_ppm;
		const dreal ifeq_ppz = no1/ffeq_ppz;
		const dreal ifeq_ppp = no1/ffeq_ppp;
		
		//raw moments used in tensors
		const dreal M200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal M020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal M002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];

		const dreal M110 = KS.f[mmm] + KS.f[mmz] + KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] + KS.f[ppp];
		const dreal M101 = KS.f[mmm] - KS.f[mmp] + KS.f[mzm] - KS.f[mzp] + KS.f[mpm] - KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] + KS.f[ppp];
		const dreal M011 = KS.f[mmm] - KS.f[mmp] - KS.f[mpm] + KS.f[mpp] + KS.f[zmm] - KS.f[zmp] - KS.f[zpm] + KS.f[zpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] + KS.f[ppp];
		const dreal M111 = KS.f[ppp] + KS.f[mmp] + KS.f[mpm] - KS.f[mpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] - KS.f[mmm];

		const dreal M201 = KS.f[ppp] + KS.f[mmp] - KS.f[mzm] + KS.f[mzp] - KS.f[mpm] + KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] - KS.f[mmm];
		const dreal M102 = KS.f[ppp] - KS.f[mmp] - KS.f[mzm] - KS.f[mzp] - KS.f[mpm] - KS.f[mpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] - KS.f[mmm];
		const dreal M210 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] + KS.f[mpm] + KS.f[mpz] + KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M120 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] + KS.f[pmm] + KS.f[pmz] + KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M021 = KS.f[ppp] + KS.f[mmp] - KS.f[mpm] + KS.f[mpp] - KS.f[zmm] + KS.f[zmp] - KS.f[zpm] + KS.f[zpp] - KS.f[pmm] + KS.f[pmp] - KS.f[ppm] - KS.f[mmm];
		const dreal M012 = KS.f[ppp] - KS.f[mmp] + KS.f[mpm] + KS.f[mpp] - KS.f[zmm] - KS.f[zmp] + KS.f[zpm] + KS.f[zpp] - KS.f[pmm] - KS.f[pmp] + KS.f[ppm] - KS.f[mmm];

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D
		const dreal	Dsmmm = Dmmm - eDmmm;
		const dreal Dsmmz = Dmmz - eDmmz;
		const dreal Dsmmp = Dmmp - eDmmp;
		const dreal Dsmzm = Dmzm - eDmzm;
		const dreal Dsmzz = Dmzz - eDmzz;
		const dreal Dsmzp = Dmzp - eDmzp;
		const dreal Dsmpm = Dmpm - eDmpm;
		const dreal Dsmpz = Dmpz - eDmpz;
		const dreal Dsmpp = Dmpp - eDmpp;
		const dreal Dszmm = Dzmm - eDzmm;
		const dreal Dszmz = Dzmz - eDzmz;
		const dreal Dszmp = Dzmp - eDzmp;
		const dreal Dszzm = Dzzm - eDzzm;
		const dreal Dszzz = Dzzz - eDzzz;
		const dreal Dszzp = Dzzp - eDzzp;
		const dreal Dszpm = Dzpm - eDzpm;
		const dreal Dszpz = Dzpz - eDzpz;
		const dreal Dszpp = Dzpp - eDzpp;
		const dreal Dspmm = Dpmm - eDpmm;
		const dreal Dspmz = Dpmz - eDpmz;
		const dreal Dspmp = Dpmp - eDpmp;
		const dreal Dspzm = Dpzm - eDpzm;
		const dreal Dspzz = Dpzz - eDpzz;
		const dreal Dspzp = Dpzp - eDpzp;
		const dreal Dsppm = Dppm - eDppm;
		const dreal Dsppz = Dppz - eDppz;
		const dreal Dsppp = Dppp - eDppp;

		// delta of high order part: delta h_i = h_i - h_i^eq = f_i - f_i^eq - delta s_i
		const dreal	Dhmmm = KS.f[mmm] - ffeq_mmm - Dsmmm;
		const dreal Dhmmz = KS.f[mmz] - ffeq_mmz - Dsmmz;
		const dreal Dhmmp = KS.f[mmp] - ffeq_mmp - Dsmmp;
		const dreal Dhmzm = KS.f[mzm] - ffeq_mzm - Dsmzm;
		const dreal Dhmzz = KS.f[mzz] - ffeq_mzz - Dsmzz;
		const dreal Dhmzp = KS.f[mzp] - ffeq_mzp - Dsmzp;
		const dreal Dhmpm = KS.f[mpm] - ffeq_mpm - Dsmpm;
		const dreal Dhmpz = KS.f[mpz] - ffeq_mpz - Dsmpz;
		const dreal Dhmpp = KS.f[mpp] - ffeq_mpp - Dsmpp;
		const dreal Dhzmm = KS.f[zmm] - ffeq_zmm - Dszmm;
		const dreal Dhzmz = KS.f[zmz] - ffeq_zmz - Dszmz;
		const dreal Dhzmp = KS.f[zmp] - ffeq_zmp - Dszmp;
		const dreal Dhzzm = KS.f[zzm] - ffeq_zzm - Dszzm;
		const dreal Dhzzz = KS.f[zzz] - ffeq_zzz - Dszzz;
		const dreal Dhzzp = KS.f[zzp] - ffeq_zzp - Dszzp;
		const dreal Dhzpm = KS.f[zpm] - ffeq_zpm - Dszpm;
		const dreal Dhzpz = KS.f[zpz] - ffeq_zpz - Dszpz;
		const dreal Dhzpp = KS.f[zpp] - ffeq_zpp - Dszpp;
		const dreal Dhpmm = KS.f[pmm] - ffeq_pmm - Dspmm;
		const dreal Dhpmz = KS.f[pmz] - ffeq_pmz - Dspmz;
		const dreal Dhpmp = KS.f[pmp] - ffeq_pmp - Dspmp;
		const dreal Dhpzm = KS.f[pzm] - ffeq_pzm - Dspzm;
		const dreal Dhpzz = KS.f[pzz] - ffeq_pzz - Dspzz;
		const dreal Dhpzp = KS.f[pzp] - ffeq_pzp - Dspzp;
		const dreal Dhppm = KS.f[ppm] - ffeq_ppm - Dsppm;
		const dreal Dhppz = KS.f[ppz] - ffeq_ppz - Dsppz;
		const dreal Dhppp = KS.f[ppp] - ffeq_ppp - Dsppp;

		const dreal beta = _beta;
		const dreal gamma = _gamma;

		// Forcing - BGK-like
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

		KS.f[mmm] -= beta*(no2*Dsmmm + gamma*Dhmmm) - (no1 - beta)*Smmm*ffeq_mmm;
		KS.f[mmz] -= beta*(no2*Dsmmz + gamma*Dhmmz) - (no1 - beta)*Smmz*ffeq_mmz;
		KS.f[mmp] -= beta*(no2*Dsmmp + gamma*Dhmmp) - (no1 - beta)*Smmp*ffeq_mmp;
		KS.f[mzm] -= beta*(no2*Dsmzm + gamma*Dhmzm) - (no1 - beta)*Smzm*ffeq_mzm;
		KS.f[mzz] -= beta*(no2*Dsmzz + gamma*Dhmzz) - (no1 - beta)*Smzz*ffeq_mzz;
		KS.f[mzp] -= beta*(no2*Dsmzp + gamma*Dhmzp) - (no1 - beta)*Smzp*ffeq_mzp;
		KS.f[mpm] -= beta*(no2*Dsmpm + gamma*Dhmpm) - (no1 - beta)*Smpm*ffeq_mpm;
		KS.f[mpz] -= beta*(no2*Dsmpz + gamma*Dhmpz) - (no1 - beta)*Smpz*ffeq_mpz;
		KS.f[mpp] -= beta*(no2*Dsmpp + gamma*Dhmpp) - (no1 - beta)*Smpp*ffeq_mpp;
		KS.f[zmm] -= beta*(no2*Dszmm + gamma*Dhzmm) - (no1 - beta)*Szmm*ffeq_zmm;
		KS.f[zmz] -= beta*(no2*Dszmz + gamma*Dhzmz) - (no1 - beta)*Szmz*ffeq_zmz;
		KS.f[zmp] -= beta*(no2*Dszmp + gamma*Dhzmp) - (no1 - beta)*Szmp*ffeq_zmp;
		KS.f[zzm] -= beta*(no2*Dszzm + gamma*Dhzzm) - (no1 - beta)*Szzm*ffeq_zzm;
		KS.f[zzz] -= beta*(no2*Dszzz + gamma*Dhzzz) - (no1 - beta)*Szzz*ffeq_zzz;
		KS.f[zzp] -= beta*(no2*Dszzp + gamma*Dhzzp) - (no1 - beta)*Szzp*ffeq_zzp;
		KS.f[zpm] -= beta*(no2*Dszpm + gamma*Dhzpm) - (no1 - beta)*Szpm*ffeq_zpm;
		KS.f[zpz] -= beta*(no2*Dszpz + gamma*Dhzpz) - (no1 - beta)*Szpz*ffeq_zpz;
		KS.f[zpp] -= beta*(no2*Dszpp + gamma*Dhzpp) - (no1 - beta)*Szpp*ffeq_zpp;
		KS.f[pmm] -= beta*(no2*Dspmm + gamma*Dhpmm) - (no1 - beta)*Spmm*ffeq_pmm;
		KS.f[pmz] -= beta*(no2*Dspmz + gamma*Dhpmz) - (no1 - beta)*Spmz*ffeq_pmz;
		KS.f[pmp] -= beta*(no2*Dspmp + gamma*Dhpmp) - (no1 - beta)*Spmp*ffeq_pmp;
		KS.f[pzm] -= beta*(no2*Dspzm + gamma*Dhpzm) - (no1 - beta)*Spzm*ffeq_pzm;
		KS.f[pzz] -= beta*(no2*Dspzz + gamma*Dhpzz) - (no1 - beta)*Spzz*ffeq_pzz;
		KS.f[pzp] -= beta*(no2*Dspzp + gamma*Dhpzp) - (no1 - beta)*Spzp*ffeq_pzp;
		KS.f[ppm] -= beta*(no2*Dsppm + gamma*Dhppm) - (no1 - beta)*Sppm*ffeq_ppm;
		KS.f[ppz] -= beta*(no2*Dsppz + gamma*Dhppz) - (no1 - beta)*Sppz*ffeq_ppz;
		KS.f[ppp] -= beta*(no2*Dsppp + gamma*Dhppp) - (no1 - beta)*Sppp*ffeq_ppp;
	}
};


template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_ENTROPIC<TRAITS>
>
struct LBM_KBC_N2 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_N2";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
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
		const dreal ffeq_mmm = -KS.rho*Xm*Ym*Zm;
		const dreal ffeq_mmz = -KS.rho*Xm*Ym*Zz;
		const dreal ffeq_mmp = -KS.rho*Xm*Ym*Zp;
		const dreal ffeq_mzm = -KS.rho*Xm*Yz*Zm;
		const dreal ffeq_mzz = -KS.rho*Xm*Yz*Zz;
		const dreal ffeq_mzp = -KS.rho*Xm*Yz*Zp;
		const dreal ffeq_mpm = -KS.rho*Xm*Yp*Zm;
		const dreal ffeq_mpz = -KS.rho*Xm*Yp*Zz;
		const dreal ffeq_mpp = -KS.rho*Xm*Yp*Zp;
		const dreal ffeq_zmm = -KS.rho*Xz*Ym*Zm;
		const dreal ffeq_zmz = -KS.rho*Xz*Ym*Zz;
		const dreal ffeq_zmp = -KS.rho*Xz*Ym*Zp;
		const dreal ffeq_zzm = -KS.rho*Xz*Yz*Zm;
		const dreal ffeq_zzz = -KS.rho*Xz*Yz*Zz;
		const dreal ffeq_zzp = -KS.rho*Xz*Yz*Zp;
		const dreal ffeq_zpm = -KS.rho*Xz*Yp*Zm;
		const dreal ffeq_zpz = -KS.rho*Xz*Yp*Zz;
		const dreal ffeq_zpp = -KS.rho*Xz*Yp*Zp;
		const dreal ffeq_pmm = -KS.rho*Xp*Ym*Zm;
		const dreal ffeq_pmz = -KS.rho*Xp*Ym*Zz;
		const dreal ffeq_pmp = -KS.rho*Xp*Ym*Zp;
		const dreal ffeq_pzm = -KS.rho*Xp*Yz*Zm;
		const dreal ffeq_pzz = -KS.rho*Xp*Yz*Zz;
		const dreal ffeq_pzp = -KS.rho*Xp*Yz*Zp;
		const dreal ffeq_ppm = -KS.rho*Xp*Yp*Zm;
		const dreal ffeq_ppz = -KS.rho*Xp*Yp*Zz;
		const dreal ffeq_ppp = -KS.rho*Xp*Yp*Zp;

		const dreal	ifeq_mmm = no1/ffeq_mmm;
		const dreal ifeq_mmz = no1/ffeq_mmz;
		const dreal ifeq_mmp = no1/ffeq_mmp;
		const dreal ifeq_mzm = no1/ffeq_mzm;
		const dreal ifeq_mzz = no1/ffeq_mzz;
		const dreal ifeq_mzp = no1/ffeq_mzp;
		const dreal ifeq_mpm = no1/ffeq_mpm;
		const dreal ifeq_mpz = no1/ffeq_mpz;
		const dreal ifeq_mpp = no1/ffeq_mpp;
		const dreal ifeq_zmm = no1/ffeq_zmm;
		const dreal ifeq_zmz = no1/ffeq_zmz;
		const dreal ifeq_zmp = no1/ffeq_zmp;
		const dreal ifeq_zzm = no1/ffeq_zzm;
		const dreal ifeq_zzz = no1/ffeq_zzz;
		const dreal ifeq_zzp = no1/ffeq_zzp;
		const dreal ifeq_zpm = no1/ffeq_zpm;
		const dreal ifeq_zpz = no1/ffeq_zpz;
		const dreal ifeq_zpp = no1/ffeq_zpp;
		const dreal ifeq_pmm = no1/ffeq_pmm;
		const dreal ifeq_pmz = no1/ffeq_pmz;
		const dreal ifeq_pmp = no1/ffeq_pmp;
		const dreal ifeq_pzm = no1/ffeq_pzm;
		const dreal ifeq_pzz = no1/ffeq_pzz;
		const dreal ifeq_pzp = no1/ffeq_pzp;
		const dreal ifeq_ppm = no1/ffeq_ppm;
		const dreal ifeq_ppz = no1/ffeq_ppz;
		const dreal ifeq_ppp = no1/ffeq_ppp;
		
		//raw moments used in tensors
		const dreal M200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal M020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal M002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];

		const dreal M110 = KS.f[mmm] + KS.f[mmz] + KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] + KS.f[ppp];
		const dreal M101 = KS.f[mmm] - KS.f[mmp] + KS.f[mzm] - KS.f[mzp] + KS.f[mpm] - KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] + KS.f[ppp];
		const dreal M011 = KS.f[mmm] - KS.f[mmp] - KS.f[mpm] + KS.f[mpp] + KS.f[zmm] - KS.f[zmp] - KS.f[zpm] + KS.f[zpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] + KS.f[ppp];
		const dreal M111 = KS.f[ppp] + KS.f[mmp] + KS.f[mpm] - KS.f[mpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] - KS.f[mmm];

		const dreal M201 = KS.f[ppp] + KS.f[mmp] - KS.f[mzm] + KS.f[mzp] - KS.f[mpm] + KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] - KS.f[mmm];
		const dreal M102 = KS.f[ppp] - KS.f[mmp] - KS.f[mzm] - KS.f[mzp] - KS.f[mpm] - KS.f[mpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] - KS.f[mmm];
		const dreal M210 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] + KS.f[mpm] + KS.f[mpz] + KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M120 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] + KS.f[pmm] + KS.f[pmz] + KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M021 = KS.f[ppp] + KS.f[mmp] - KS.f[mpm] + KS.f[mpp] - KS.f[zmm] + KS.f[zmp] - KS.f[zpm] + KS.f[zpp] - KS.f[pmm] + KS.f[pmp] - KS.f[ppm] - KS.f[mmm];
		const dreal M012 = KS.f[ppp] - KS.f[mmp] + KS.f[mpm] + KS.f[mpp] - KS.f[zmm] - KS.f[zmp] + KS.f[zpm] + KS.f[zpp] - KS.f[pmm] - KS.f[pmp] + KS.f[ppm] - KS.f[mmm];

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D + T
		const dreal	Dsmmm = Dmmm - eDmmm + Tmmm - eTmmm;
		const dreal Dsmmz = Dmmz - eDmmz + Tmmz - eTmmz;
		const dreal Dsmmp = Dmmp - eDmmp + Tmmp - eTmmp;
		const dreal Dsmzm = Dmzm - eDmzm + Tmzm - eTmzm;
		const dreal Dsmzz = Dmzz - eDmzz + Tmzz - eTmzz;
		const dreal Dsmzp = Dmzp - eDmzp + Tmzp - eTmzp;
		const dreal Dsmpm = Dmpm - eDmpm + Tmpm - eTmpm;
		const dreal Dsmpz = Dmpz - eDmpz + Tmpz - eTmpz;
		const dreal Dsmpp = Dmpp - eDmpp + Tmpp - eTmpp;
		const dreal Dszmm = Dzmm - eDzmm + Tzmm - eTzmm;
		const dreal Dszmz = Dzmz - eDzmz + Tzmz - eTzmz;
		const dreal Dszmp = Dzmp - eDzmp + Tzmp - eTzmp;
		const dreal Dszzm = Dzzm - eDzzm + Tzzm - eTzzm;
		const dreal Dszzz = Dzzz - eDzzz + Tzzz - eTzzz;
		const dreal Dszzp = Dzzp - eDzzp + Tzzp - eTzzp;
		const dreal Dszpm = Dzpm - eDzpm + Tzpm - eTzpm;
		const dreal Dszpz = Dzpz - eDzpz + Tzpz - eTzpz;
		const dreal Dszpp = Dzpp - eDzpp + Tzpp - eTzpp;
		const dreal Dspmm = Dpmm - eDpmm + Tpmm - eTpmm;
		const dreal Dspmz = Dpmz - eDpmz + Tpmz - eTpmz;
		const dreal Dspmp = Dpmp - eDpmp + Tpmp - eTpmp;
		const dreal Dspzm = Dpzm - eDpzm + Tpzm - eTpzm;
		const dreal Dspzz = Dpzz - eDpzz + Tpzz - eTpzz;
		const dreal Dspzp = Dpzp - eDpzp + Tpzp - eTpzp;
		const dreal Dsppm = Dppm - eDppm + Tppm - eTppm;
		const dreal Dsppz = Dppz - eDppz + Tppz - eTppz;
		const dreal Dsppp = Dppp - eDppp + Tppp - eTppp;

		// delta of high order part: delta h_i = h_i - h_i^eq = f_i - f_i^eq - delta s_i
		const dreal	Dhmmm = KS.f[mmm] - ffeq_mmm - Dsmmm;
		const dreal Dhmmz = KS.f[mmz] - ffeq_mmz - Dsmmz;
		const dreal Dhmmp = KS.f[mmp] - ffeq_mmp - Dsmmp;
		const dreal Dhmzm = KS.f[mzm] - ffeq_mzm - Dsmzm;
		const dreal Dhmzz = KS.f[mzz] - ffeq_mzz - Dsmzz;
		const dreal Dhmzp = KS.f[mzp] - ffeq_mzp - Dsmzp;
		const dreal Dhmpm = KS.f[mpm] - ffeq_mpm - Dsmpm;
		const dreal Dhmpz = KS.f[mpz] - ffeq_mpz - Dsmpz;
		const dreal Dhmpp = KS.f[mpp] - ffeq_mpp - Dsmpp;
		const dreal Dhzmm = KS.f[zmm] - ffeq_zmm - Dszmm;
		const dreal Dhzmz = KS.f[zmz] - ffeq_zmz - Dszmz;
		const dreal Dhzmp = KS.f[zmp] - ffeq_zmp - Dszmp;
		const dreal Dhzzm = KS.f[zzm] - ffeq_zzm - Dszzm;
		const dreal Dhzzz = KS.f[zzz] - ffeq_zzz - Dszzz;
		const dreal Dhzzp = KS.f[zzp] - ffeq_zzp - Dszzp;
		const dreal Dhzpm = KS.f[zpm] - ffeq_zpm - Dszpm;
		const dreal Dhzpz = KS.f[zpz] - ffeq_zpz - Dszpz;
		const dreal Dhzpp = KS.f[zpp] - ffeq_zpp - Dszpp;
		const dreal Dhpmm = KS.f[pmm] - ffeq_pmm - Dspmm;
		const dreal Dhpmz = KS.f[pmz] - ffeq_pmz - Dspmz;
		const dreal Dhpmp = KS.f[pmp] - ffeq_pmp - Dspmp;
		const dreal Dhpzm = KS.f[pzm] - ffeq_pzm - Dspzm;
		const dreal Dhpzz = KS.f[pzz] - ffeq_pzz - Dspzz;
		const dreal Dhpzp = KS.f[pzp] - ffeq_pzp - Dspzp;
		const dreal Dhppm = KS.f[ppm] - ffeq_ppm - Dsppm;
		const dreal Dhppz = KS.f[ppz] - ffeq_ppz - Dsppz;
		const dreal Dhppp = KS.f[ppp] - ffeq_ppp - Dsppp;

		const dreal beta = _beta;
		const dreal gamma = _gamma;

		// Forcing - BGK-like
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

		KS.f[mmm] -= beta*(no2*Dsmmm + gamma*Dhmmm) - (no1 - beta)*Smmm*ffeq_mmm;
		KS.f[mmz] -= beta*(no2*Dsmmz + gamma*Dhmmz) - (no1 - beta)*Smmz*ffeq_mmz;
		KS.f[mmp] -= beta*(no2*Dsmmp + gamma*Dhmmp) - (no1 - beta)*Smmp*ffeq_mmp;
		KS.f[mzm] -= beta*(no2*Dsmzm + gamma*Dhmzm) - (no1 - beta)*Smzm*ffeq_mzm;
		KS.f[mzz] -= beta*(no2*Dsmzz + gamma*Dhmzz) - (no1 - beta)*Smzz*ffeq_mzz;
		KS.f[mzp] -= beta*(no2*Dsmzp + gamma*Dhmzp) - (no1 - beta)*Smzp*ffeq_mzp;
		KS.f[mpm] -= beta*(no2*Dsmpm + gamma*Dhmpm) - (no1 - beta)*Smpm*ffeq_mpm;
		KS.f[mpz] -= beta*(no2*Dsmpz + gamma*Dhmpz) - (no1 - beta)*Smpz*ffeq_mpz;
		KS.f[mpp] -= beta*(no2*Dsmpp + gamma*Dhmpp) - (no1 - beta)*Smpp*ffeq_mpp;
		KS.f[zmm] -= beta*(no2*Dszmm + gamma*Dhzmm) - (no1 - beta)*Szmm*ffeq_zmm;
		KS.f[zmz] -= beta*(no2*Dszmz + gamma*Dhzmz) - (no1 - beta)*Szmz*ffeq_zmz;
		KS.f[zmp] -= beta*(no2*Dszmp + gamma*Dhzmp) - (no1 - beta)*Szmp*ffeq_zmp;
		KS.f[zzm] -= beta*(no2*Dszzm + gamma*Dhzzm) - (no1 - beta)*Szzm*ffeq_zzm;
		KS.f[zzz] -= beta*(no2*Dszzz + gamma*Dhzzz) - (no1 - beta)*Szzz*ffeq_zzz;
		KS.f[zzp] -= beta*(no2*Dszzp + gamma*Dhzzp) - (no1 - beta)*Szzp*ffeq_zzp;
		KS.f[zpm] -= beta*(no2*Dszpm + gamma*Dhzpm) - (no1 - beta)*Szpm*ffeq_zpm;
		KS.f[zpz] -= beta*(no2*Dszpz + gamma*Dhzpz) - (no1 - beta)*Szpz*ffeq_zpz;
		KS.f[zpp] -= beta*(no2*Dszpp + gamma*Dhzpp) - (no1 - beta)*Szpp*ffeq_zpp;
		KS.f[pmm] -= beta*(no2*Dspmm + gamma*Dhpmm) - (no1 - beta)*Spmm*ffeq_pmm;
		KS.f[pmz] -= beta*(no2*Dspmz + gamma*Dhpmz) - (no1 - beta)*Spmz*ffeq_pmz;
		KS.f[pmp] -= beta*(no2*Dspmp + gamma*Dhpmp) - (no1 - beta)*Spmp*ffeq_pmp;
		KS.f[pzm] -= beta*(no2*Dspzm + gamma*Dhpzm) - (no1 - beta)*Spzm*ffeq_pzm;
		KS.f[pzz] -= beta*(no2*Dspzz + gamma*Dhpzz) - (no1 - beta)*Spzz*ffeq_pzz;
		KS.f[pzp] -= beta*(no2*Dspzp + gamma*Dhpzp) - (no1 - beta)*Spzp*ffeq_pzp;
		KS.f[ppm] -= beta*(no2*Dsppm + gamma*Dhppm) - (no1 - beta)*Sppm*ffeq_ppm;
		KS.f[ppz] -= beta*(no2*Dsppz + gamma*Dhppz) - (no1 - beta)*Sppz*ffeq_ppz;
		KS.f[ppp] -= beta*(no2*Dsppp + gamma*Dhppp) - (no1 - beta)*Sppp*ffeq_ppp;
	}
};


template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_ENTROPIC<TRAITS>
>
struct LBM_KBC_N3 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_N3";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
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
		const dreal ffeq_mmm = -KS.rho*Xm*Ym*Zm;
		const dreal ffeq_mmz = -KS.rho*Xm*Ym*Zz;
		const dreal ffeq_mmp = -KS.rho*Xm*Ym*Zp;
		const dreal ffeq_mzm = -KS.rho*Xm*Yz*Zm;
		const dreal ffeq_mzz = -KS.rho*Xm*Yz*Zz;
		const dreal ffeq_mzp = -KS.rho*Xm*Yz*Zp;
		const dreal ffeq_mpm = -KS.rho*Xm*Yp*Zm;
		const dreal ffeq_mpz = -KS.rho*Xm*Yp*Zz;
		const dreal ffeq_mpp = -KS.rho*Xm*Yp*Zp;
		const dreal ffeq_zmm = -KS.rho*Xz*Ym*Zm;
		const dreal ffeq_zmz = -KS.rho*Xz*Ym*Zz;
		const dreal ffeq_zmp = -KS.rho*Xz*Ym*Zp;
		const dreal ffeq_zzm = -KS.rho*Xz*Yz*Zm;
		const dreal ffeq_zzz = -KS.rho*Xz*Yz*Zz;
		const dreal ffeq_zzp = -KS.rho*Xz*Yz*Zp;
		const dreal ffeq_zpm = -KS.rho*Xz*Yp*Zm;
		const dreal ffeq_zpz = -KS.rho*Xz*Yp*Zz;
		const dreal ffeq_zpp = -KS.rho*Xz*Yp*Zp;
		const dreal ffeq_pmm = -KS.rho*Xp*Ym*Zm;
		const dreal ffeq_pmz = -KS.rho*Xp*Ym*Zz;
		const dreal ffeq_pmp = -KS.rho*Xp*Ym*Zp;
		const dreal ffeq_pzm = -KS.rho*Xp*Yz*Zm;
		const dreal ffeq_pzz = -KS.rho*Xp*Yz*Zz;
		const dreal ffeq_pzp = -KS.rho*Xp*Yz*Zp;
		const dreal ffeq_ppm = -KS.rho*Xp*Yp*Zm;
		const dreal ffeq_ppz = -KS.rho*Xp*Yp*Zz;
		const dreal ffeq_ppp = -KS.rho*Xp*Yp*Zp;

		const dreal	ifeq_mmm = no1/ffeq_mmm;
		const dreal ifeq_mmz = no1/ffeq_mmz;
		const dreal ifeq_mmp = no1/ffeq_mmp;
		const dreal ifeq_mzm = no1/ffeq_mzm;
		const dreal ifeq_mzz = no1/ffeq_mzz;
		const dreal ifeq_mzp = no1/ffeq_mzp;
		const dreal ifeq_mpm = no1/ffeq_mpm;
		const dreal ifeq_mpz = no1/ffeq_mpz;
		const dreal ifeq_mpp = no1/ffeq_mpp;
		const dreal ifeq_zmm = no1/ffeq_zmm;
		const dreal ifeq_zmz = no1/ffeq_zmz;
		const dreal ifeq_zmp = no1/ffeq_zmp;
		const dreal ifeq_zzm = no1/ffeq_zzm;
		const dreal ifeq_zzz = no1/ffeq_zzz;
		const dreal ifeq_zzp = no1/ffeq_zzp;
		const dreal ifeq_zpm = no1/ffeq_zpm;
		const dreal ifeq_zpz = no1/ffeq_zpz;
		const dreal ifeq_zpp = no1/ffeq_zpp;
		const dreal ifeq_pmm = no1/ffeq_pmm;
		const dreal ifeq_pmz = no1/ffeq_pmz;
		const dreal ifeq_pmp = no1/ffeq_pmp;
		const dreal ifeq_pzm = no1/ffeq_pzm;
		const dreal ifeq_pzz = no1/ffeq_pzz;
		const dreal ifeq_pzp = no1/ffeq_pzp;
		const dreal ifeq_ppm = no1/ffeq_ppm;
		const dreal ifeq_ppz = no1/ffeq_ppz;
		const dreal ifeq_ppp = no1/ffeq_ppp;
		
		//raw moments used in tensors
		const dreal M200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal M020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal M002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];

		const dreal M110 = KS.f[mmm] + KS.f[mmz] + KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] + KS.f[ppp];
		const dreal M101 = KS.f[mmm] - KS.f[mmp] + KS.f[mzm] - KS.f[mzp] + KS.f[mpm] - KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] + KS.f[ppp];
		const dreal M011 = KS.f[mmm] - KS.f[mmp] - KS.f[mpm] + KS.f[mpp] + KS.f[zmm] - KS.f[zmp] - KS.f[zpm] + KS.f[zpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] + KS.f[ppp];
		const dreal M111 = KS.f[ppp] + KS.f[mmp] + KS.f[mpm] - KS.f[mpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] - KS.f[mmm];

		const dreal M201 = KS.f[ppp] + KS.f[mmp] - KS.f[mzm] + KS.f[mzp] - KS.f[mpm] + KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] - KS.f[mmm];
		const dreal M102 = KS.f[ppp] - KS.f[mmp] - KS.f[mzm] - KS.f[mzp] - KS.f[mpm] - KS.f[mpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] - KS.f[mmm];
		const dreal M210 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] + KS.f[mpm] + KS.f[mpz] + KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M120 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] + KS.f[pmm] + KS.f[pmz] + KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M021 = KS.f[ppp] + KS.f[mmp] - KS.f[mpm] + KS.f[mpp] - KS.f[zmm] + KS.f[zmp] - KS.f[zpm] + KS.f[zpp] - KS.f[pmm] + KS.f[pmp] - KS.f[ppm] - KS.f[mmm];
		const dreal M012 = KS.f[ppp] - KS.f[mmp] + KS.f[mpm] + KS.f[mpp] - KS.f[zmm] - KS.f[zmp] + KS.f[zpm] + KS.f[zpp] - KS.f[pmm] - KS.f[pmp] + KS.f[ppm] - KS.f[mmm];

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D + Q
		const dreal	Dsmmm = Dmmm - eDmmm + Qmmm - eQmmm;
		const dreal Dsmmz = Dmmz - eDmmz + Qmmz - eQmmz;
		const dreal Dsmmp = Dmmp - eDmmp + Qmmp - eQmmp;
		const dreal Dsmzm = Dmzm - eDmzm + Qmzm - eQmzm;
		const dreal Dsmzz = Dmzz - eDmzz + Qmzz - eQmzz;
		const dreal Dsmzp = Dmzp - eDmzp + Qmzp - eQmzp;
		const dreal Dsmpm = Dmpm - eDmpm + Qmpm - eQmpm;
		const dreal Dsmpz = Dmpz - eDmpz + Qmpz - eQmpz;
		const dreal Dsmpp = Dmpp - eDmpp + Qmpp - eQmpp;
		const dreal Dszmm = Dzmm - eDzmm + Qzmm - eQzmm;
		const dreal Dszmz = Dzmz - eDzmz + Qzmz - eQzmz;
		const dreal Dszmp = Dzmp - eDzmp + Qzmp - eQzmp;
		const dreal Dszzm = Dzzm - eDzzm + Qzzm - eQzzm;
		const dreal Dszzz = Dzzz - eDzzz + Qzzz - eQzzz;
		const dreal Dszzp = Dzzp - eDzzp + Qzzp - eQzzp;
		const dreal Dszpm = Dzpm - eDzpm + Qzpm - eQzpm;
		const dreal Dszpz = Dzpz - eDzpz + Qzpz - eQzpz;
		const dreal Dszpp = Dzpp - eDzpp + Qzpp - eQzpp;
		const dreal Dspmm = Dpmm - eDpmm + Qpmm - eQpmm;
		const dreal Dspmz = Dpmz - eDpmz + Qpmz - eQpmz;
		const dreal Dspmp = Dpmp - eDpmp + Qpmp - eQpmp;
		const dreal Dspzm = Dpzm - eDpzm + Qpzm - eQpzm;
		const dreal Dspzz = Dpzz - eDpzz + Qpzz - eQpzz;
		const dreal Dspzp = Dpzp - eDpzp + Qpzp - eQpzp;
		const dreal Dsppm = Dppm - eDppm + Qppm - eQppm;
		const dreal Dsppz = Dppz - eDppz + Qppz - eQppz;
		const dreal Dsppp = Dppp - eDppp + Qppp - eQppp;

		// delta of high order part: delta h_i = h_i - h_i^eq = f_i - f_i^eq - delta s_i
		const dreal	Dhmmm = KS.f[mmm] - ffeq_mmm - Dsmmm;
		const dreal Dhmmz = KS.f[mmz] - ffeq_mmz - Dsmmz;
		const dreal Dhmmp = KS.f[mmp] - ffeq_mmp - Dsmmp;
		const dreal Dhmzm = KS.f[mzm] - ffeq_mzm - Dsmzm;
		const dreal Dhmzz = KS.f[mzz] - ffeq_mzz - Dsmzz;
		const dreal Dhmzp = KS.f[mzp] - ffeq_mzp - Dsmzp;
		const dreal Dhmpm = KS.f[mpm] - ffeq_mpm - Dsmpm;
		const dreal Dhmpz = KS.f[mpz] - ffeq_mpz - Dsmpz;
		const dreal Dhmpp = KS.f[mpp] - ffeq_mpp - Dsmpp;
		const dreal Dhzmm = KS.f[zmm] - ffeq_zmm - Dszmm;
		const dreal Dhzmz = KS.f[zmz] - ffeq_zmz - Dszmz;
		const dreal Dhzmp = KS.f[zmp] - ffeq_zmp - Dszmp;
		const dreal Dhzzm = KS.f[zzm] - ffeq_zzm - Dszzm;
		const dreal Dhzzz = KS.f[zzz] - ffeq_zzz - Dszzz;
		const dreal Dhzzp = KS.f[zzp] - ffeq_zzp - Dszzp;
		const dreal Dhzpm = KS.f[zpm] - ffeq_zpm - Dszpm;
		const dreal Dhzpz = KS.f[zpz] - ffeq_zpz - Dszpz;
		const dreal Dhzpp = KS.f[zpp] - ffeq_zpp - Dszpp;
		const dreal Dhpmm = KS.f[pmm] - ffeq_pmm - Dspmm;
		const dreal Dhpmz = KS.f[pmz] - ffeq_pmz - Dspmz;
		const dreal Dhpmp = KS.f[pmp] - ffeq_pmp - Dspmp;
		const dreal Dhpzm = KS.f[pzm] - ffeq_pzm - Dspzm;
		const dreal Dhpzz = KS.f[pzz] - ffeq_pzz - Dspzz;
		const dreal Dhpzp = KS.f[pzp] - ffeq_pzp - Dspzp;
		const dreal Dhppm = KS.f[ppm] - ffeq_ppm - Dsppm;
		const dreal Dhppz = KS.f[ppz] - ffeq_ppz - Dsppz;
		const dreal Dhppp = KS.f[ppp] - ffeq_ppp - Dsppp;

		const dreal beta = _beta;
		const dreal gamma = _gamma;

		// Forcing - BGK-like
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

		KS.f[mmm] -= beta*(no2*Dsmmm + gamma*Dhmmm) - (no1 - beta)*Smmm*ffeq_mmm;
		KS.f[mmz] -= beta*(no2*Dsmmz + gamma*Dhmmz) - (no1 - beta)*Smmz*ffeq_mmz;
		KS.f[mmp] -= beta*(no2*Dsmmp + gamma*Dhmmp) - (no1 - beta)*Smmp*ffeq_mmp;
		KS.f[mzm] -= beta*(no2*Dsmzm + gamma*Dhmzm) - (no1 - beta)*Smzm*ffeq_mzm;
		KS.f[mzz] -= beta*(no2*Dsmzz + gamma*Dhmzz) - (no1 - beta)*Smzz*ffeq_mzz;
		KS.f[mzp] -= beta*(no2*Dsmzp + gamma*Dhmzp) - (no1 - beta)*Smzp*ffeq_mzp;
		KS.f[mpm] -= beta*(no2*Dsmpm + gamma*Dhmpm) - (no1 - beta)*Smpm*ffeq_mpm;
		KS.f[mpz] -= beta*(no2*Dsmpz + gamma*Dhmpz) - (no1 - beta)*Smpz*ffeq_mpz;
		KS.f[mpp] -= beta*(no2*Dsmpp + gamma*Dhmpp) - (no1 - beta)*Smpp*ffeq_mpp;
		KS.f[zmm] -= beta*(no2*Dszmm + gamma*Dhzmm) - (no1 - beta)*Szmm*ffeq_zmm;
		KS.f[zmz] -= beta*(no2*Dszmz + gamma*Dhzmz) - (no1 - beta)*Szmz*ffeq_zmz;
		KS.f[zmp] -= beta*(no2*Dszmp + gamma*Dhzmp) - (no1 - beta)*Szmp*ffeq_zmp;
		KS.f[zzm] -= beta*(no2*Dszzm + gamma*Dhzzm) - (no1 - beta)*Szzm*ffeq_zzm;
		KS.f[zzz] -= beta*(no2*Dszzz + gamma*Dhzzz) - (no1 - beta)*Szzz*ffeq_zzz;
		KS.f[zzp] -= beta*(no2*Dszzp + gamma*Dhzzp) - (no1 - beta)*Szzp*ffeq_zzp;
		KS.f[zpm] -= beta*(no2*Dszpm + gamma*Dhzpm) - (no1 - beta)*Szpm*ffeq_zpm;
		KS.f[zpz] -= beta*(no2*Dszpz + gamma*Dhzpz) - (no1 - beta)*Szpz*ffeq_zpz;
		KS.f[zpp] -= beta*(no2*Dszpp + gamma*Dhzpp) - (no1 - beta)*Szpp*ffeq_zpp;
		KS.f[pmm] -= beta*(no2*Dspmm + gamma*Dhpmm) - (no1 - beta)*Spmm*ffeq_pmm;
		KS.f[pmz] -= beta*(no2*Dspmz + gamma*Dhpmz) - (no1 - beta)*Spmz*ffeq_pmz;
		KS.f[pmp] -= beta*(no2*Dspmp + gamma*Dhpmp) - (no1 - beta)*Spmp*ffeq_pmp;
		KS.f[pzm] -= beta*(no2*Dspzm + gamma*Dhpzm) - (no1 - beta)*Spzm*ffeq_pzm;
		KS.f[pzz] -= beta*(no2*Dspzz + gamma*Dhpzz) - (no1 - beta)*Spzz*ffeq_pzz;
		KS.f[pzp] -= beta*(no2*Dspzp + gamma*Dhpzp) - (no1 - beta)*Spzp*ffeq_pzp;
		KS.f[ppm] -= beta*(no2*Dsppm + gamma*Dhppm) - (no1 - beta)*Sppm*ffeq_ppm;
		KS.f[ppz] -= beta*(no2*Dsppz + gamma*Dhppz) - (no1 - beta)*Sppz*ffeq_ppz;
		KS.f[ppp] -= beta*(no2*Dsppp + gamma*Dhppp) - (no1 - beta)*Sppp*ffeq_ppp;
	}
};

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_ENTROPIC<TRAITS>
>
struct LBM_KBC_N4 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_N4";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
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
		const dreal ffeq_mmm = -KS.rho*Xm*Ym*Zm;
		const dreal ffeq_mmz = -KS.rho*Xm*Ym*Zz;
		const dreal ffeq_mmp = -KS.rho*Xm*Ym*Zp;
		const dreal ffeq_mzm = -KS.rho*Xm*Yz*Zm;
		const dreal ffeq_mzz = -KS.rho*Xm*Yz*Zz;
		const dreal ffeq_mzp = -KS.rho*Xm*Yz*Zp;
		const dreal ffeq_mpm = -KS.rho*Xm*Yp*Zm;
		const dreal ffeq_mpz = -KS.rho*Xm*Yp*Zz;
		const dreal ffeq_mpp = -KS.rho*Xm*Yp*Zp;
		const dreal ffeq_zmm = -KS.rho*Xz*Ym*Zm;
		const dreal ffeq_zmz = -KS.rho*Xz*Ym*Zz;
		const dreal ffeq_zmp = -KS.rho*Xz*Ym*Zp;
		const dreal ffeq_zzm = -KS.rho*Xz*Yz*Zm;
		const dreal ffeq_zzz = -KS.rho*Xz*Yz*Zz;
		const dreal ffeq_zzp = -KS.rho*Xz*Yz*Zp;
		const dreal ffeq_zpm = -KS.rho*Xz*Yp*Zm;
		const dreal ffeq_zpz = -KS.rho*Xz*Yp*Zz;
		const dreal ffeq_zpp = -KS.rho*Xz*Yp*Zp;
		const dreal ffeq_pmm = -KS.rho*Xp*Ym*Zm;
		const dreal ffeq_pmz = -KS.rho*Xp*Ym*Zz;
		const dreal ffeq_pmp = -KS.rho*Xp*Ym*Zp;
		const dreal ffeq_pzm = -KS.rho*Xp*Yz*Zm;
		const dreal ffeq_pzz = -KS.rho*Xp*Yz*Zz;
		const dreal ffeq_pzp = -KS.rho*Xp*Yz*Zp;
		const dreal ffeq_ppm = -KS.rho*Xp*Yp*Zm;
		const dreal ffeq_ppz = -KS.rho*Xp*Yp*Zz;
		const dreal ffeq_ppp = -KS.rho*Xp*Yp*Zp;

		const dreal	ifeq_mmm = no1/ffeq_mmm;
		const dreal ifeq_mmz = no1/ffeq_mmz;
		const dreal ifeq_mmp = no1/ffeq_mmp;
		const dreal ifeq_mzm = no1/ffeq_mzm;
		const dreal ifeq_mzz = no1/ffeq_mzz;
		const dreal ifeq_mzp = no1/ffeq_mzp;
		const dreal ifeq_mpm = no1/ffeq_mpm;
		const dreal ifeq_mpz = no1/ffeq_mpz;
		const dreal ifeq_mpp = no1/ffeq_mpp;
		const dreal ifeq_zmm = no1/ffeq_zmm;
		const dreal ifeq_zmz = no1/ffeq_zmz;
		const dreal ifeq_zmp = no1/ffeq_zmp;
		const dreal ifeq_zzm = no1/ffeq_zzm;
		const dreal ifeq_zzz = no1/ffeq_zzz;
		const dreal ifeq_zzp = no1/ffeq_zzp;
		const dreal ifeq_zpm = no1/ffeq_zpm;
		const dreal ifeq_zpz = no1/ffeq_zpz;
		const dreal ifeq_zpp = no1/ffeq_zpp;
		const dreal ifeq_pmm = no1/ffeq_pmm;
		const dreal ifeq_pmz = no1/ffeq_pmz;
		const dreal ifeq_pmp = no1/ffeq_pmp;
		const dreal ifeq_pzm = no1/ffeq_pzm;
		const dreal ifeq_pzz = no1/ffeq_pzz;
		const dreal ifeq_pzp = no1/ffeq_pzp;
		const dreal ifeq_ppm = no1/ffeq_ppm;
		const dreal ifeq_ppz = no1/ffeq_ppz;
		const dreal ifeq_ppp = no1/ffeq_ppp;
		
		//raw moments used in tensors
		const dreal M200 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[mzm] + KS.f[mzp] + KS.f[mzz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz] + KS.f[pzm] + KS.f[pzp] + KS.f[pzz];
		const dreal M020 = KS.f[mmm] + KS.f[mmp] + KS.f[mmz] + KS.f[mpm] + KS.f[mpp] + KS.f[mpz] + KS.f[zmm] + KS.f[zmp] + KS.f[zmz] + KS.f[zpm] + KS.f[zpp] + KS.f[zpz] + KS.f[pmm] + KS.f[pmp] + KS.f[pmz] + KS.f[ppm] + KS.f[ppp] + KS.f[ppz];
		const dreal M002 = KS.f[mmm] + KS.f[mmp] + KS.f[mzm] + KS.f[mzp] + KS.f[mpm] + KS.f[mpp] + KS.f[zmm] + KS.f[zmp] + KS.f[zzm] + KS.f[zzp] + KS.f[zpm] + KS.f[zpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] + KS.f[ppp];

		const dreal M110 = KS.f[mmm] + KS.f[mmz] + KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] + KS.f[ppp];
		const dreal M101 = KS.f[mmm] - KS.f[mmp] + KS.f[mzm] - KS.f[mzp] + KS.f[mpm] - KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] + KS.f[ppp];
		const dreal M011 = KS.f[mmm] - KS.f[mmp] - KS.f[mpm] + KS.f[mpp] + KS.f[zmm] - KS.f[zmp] - KS.f[zpm] + KS.f[zpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] + KS.f[ppp];
		const dreal M111 = KS.f[ppp] + KS.f[mmp] + KS.f[mpm] - KS.f[mpp] + KS.f[pmm] - KS.f[pmp] - KS.f[ppm] - KS.f[mmm];

		const dreal M201 = KS.f[ppp] + KS.f[mmp] - KS.f[mzm] + KS.f[mzp] - KS.f[mpm] + KS.f[mpp] - KS.f[pmm] + KS.f[pmp] - KS.f[pzm] + KS.f[pzp] - KS.f[ppm] - KS.f[mmm];
		const dreal M102 = KS.f[ppp] - KS.f[mmp] - KS.f[mzm] - KS.f[mzp] - KS.f[mpm] - KS.f[mpp] + KS.f[pmm] + KS.f[pmp] + KS.f[pzm] + KS.f[pzp] + KS.f[ppm] - KS.f[mmm];
		const dreal M210 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] + KS.f[mpm] + KS.f[mpz] + KS.f[mpp] - KS.f[pmm] - KS.f[pmz] - KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M120 = KS.f[ppp] - KS.f[mmz] - KS.f[mmp] - KS.f[mpm] - KS.f[mpz] - KS.f[mpp] + KS.f[pmm] + KS.f[pmz] + KS.f[pmp] + KS.f[ppm] + KS.f[ppz] - KS.f[mmm];
		const dreal M021 = KS.f[ppp] + KS.f[mmp] - KS.f[mpm] + KS.f[mpp] - KS.f[zmm] + KS.f[zmp] - KS.f[zpm] + KS.f[zpp] - KS.f[pmm] + KS.f[pmp] - KS.f[ppm] - KS.f[mmm];
		const dreal M012 = KS.f[ppp] - KS.f[mmp] + KS.f[mpm] + KS.f[mpp] - KS.f[zmm] - KS.f[zmp] + KS.f[zpm] + KS.f[zpp] - KS.f[pmm] - KS.f[pmp] + KS.f[ppm] - KS.f[mmm];

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D + T + Q
		const dreal	Dsmmm = Dmmm - eDmmm + Tmmm - eTmmm + Qmmm - eQmmm;
		const dreal Dsmmz = Dmmz - eDmmz + Tmmz - eTmmz + Qmmz - eQmmz;
		const dreal Dsmmp = Dmmp - eDmmp + Tmmp - eTmmp + Qmmp - eQmmp;
		const dreal Dsmzm = Dmzm - eDmzm + Tmzm - eTmzm + Qmzm - eQmzm;
		const dreal Dsmzz = Dmzz - eDmzz + Tmzz - eTmzz + Qmzz - eQmzz;
		const dreal Dsmzp = Dmzp - eDmzp + Tmzp - eTmzp + Qmzp - eQmzp;
		const dreal Dsmpm = Dmpm - eDmpm + Tmpm - eTmpm + Qmpm - eQmpm;
		const dreal Dsmpz = Dmpz - eDmpz + Tmpz - eTmpz + Qmpz - eQmpz;
		const dreal Dsmpp = Dmpp - eDmpp + Tmpp - eTmpp + Qmpp - eQmpp;
		const dreal Dszmm = Dzmm - eDzmm + Tzmm - eTzmm + Qzmm - eQzmm;
		const dreal Dszmz = Dzmz - eDzmz + Tzmz - eTzmz + Qzmz - eQzmz;
		const dreal Dszmp = Dzmp - eDzmp + Tzmp - eTzmp + Qzmp - eQzmp;
		const dreal Dszzm = Dzzm - eDzzm + Tzzm - eTzzm + Qzzm - eQzzm;
		const dreal Dszzz = Dzzz - eDzzz + Tzzz - eTzzz + Qzzz - eQzzz;
		const dreal Dszzp = Dzzp - eDzzp + Tzzp - eTzzp + Qzzp - eQzzp;
		const dreal Dszpm = Dzpm - eDzpm + Tzpm - eTzpm + Qzpm - eQzpm;
		const dreal Dszpz = Dzpz - eDzpz + Tzpz - eTzpz + Qzpz - eQzpz;
		const dreal Dszpp = Dzpp - eDzpp + Tzpp - eTzpp + Qzpp - eQzpp;
		const dreal Dspmm = Dpmm - eDpmm + Tpmm - eTpmm + Qpmm - eQpmm;
		const dreal Dspmz = Dpmz - eDpmz + Tpmz - eTpmz + Qpmz - eQpmz;
		const dreal Dspmp = Dpmp - eDpmp + Tpmp - eTpmp + Qpmp - eQpmp;
		const dreal Dspzm = Dpzm - eDpzm + Tpzm - eTpzm + Qpzm - eQpzm;
		const dreal Dspzz = Dpzz - eDpzz + Tpzz - eTpzz + Qpzz - eQpzz;
		const dreal Dspzp = Dpzp - eDpzp + Tpzp - eTpzp + Qpzp - eQpzp;
		const dreal Dsppm = Dppm - eDppm + Tppm - eTppm + Qppm - eQppm;
		const dreal Dsppz = Dppz - eDppz + Tppz - eTppz + Qppz - eQppz;
		const dreal Dsppp = Dppp - eDppp + Tppp - eTppp + Qppp - eQppp;

		// delta of high order part: delta h_i = h_i - h_i^eq = f_i - f_i^eq - delta s_i
		const dreal	Dhmmm = KS.f[mmm] - ffeq_mmm - Dsmmm;
		const dreal Dhmmz = KS.f[mmz] - ffeq_mmz - Dsmmz;
		const dreal Dhmmp = KS.f[mmp] - ffeq_mmp - Dsmmp;
		const dreal Dhmzm = KS.f[mzm] - ffeq_mzm - Dsmzm;
		const dreal Dhmzz = KS.f[mzz] - ffeq_mzz - Dsmzz;
		const dreal Dhmzp = KS.f[mzp] - ffeq_mzp - Dsmzp;
		const dreal Dhmpm = KS.f[mpm] - ffeq_mpm - Dsmpm;
		const dreal Dhmpz = KS.f[mpz] - ffeq_mpz - Dsmpz;
		const dreal Dhmpp = KS.f[mpp] - ffeq_mpp - Dsmpp;
		const dreal Dhzmm = KS.f[zmm] - ffeq_zmm - Dszmm;
		const dreal Dhzmz = KS.f[zmz] - ffeq_zmz - Dszmz;
		const dreal Dhzmp = KS.f[zmp] - ffeq_zmp - Dszmp;
		const dreal Dhzzm = KS.f[zzm] - ffeq_zzm - Dszzm;
		const dreal Dhzzz = KS.f[zzz] - ffeq_zzz - Dszzz;
		const dreal Dhzzp = KS.f[zzp] - ffeq_zzp - Dszzp;
		const dreal Dhzpm = KS.f[zpm] - ffeq_zpm - Dszpm;
		const dreal Dhzpz = KS.f[zpz] - ffeq_zpz - Dszpz;
		const dreal Dhzpp = KS.f[zpp] - ffeq_zpp - Dszpp;
		const dreal Dhpmm = KS.f[pmm] - ffeq_pmm - Dspmm;
		const dreal Dhpmz = KS.f[pmz] - ffeq_pmz - Dspmz;
		const dreal Dhpmp = KS.f[pmp] - ffeq_pmp - Dspmp;
		const dreal Dhpzm = KS.f[pzm] - ffeq_pzm - Dspzm;
		const dreal Dhpzz = KS.f[pzz] - ffeq_pzz - Dspzz;
		const dreal Dhpzp = KS.f[pzp] - ffeq_pzp - Dspzp;
		const dreal Dhppm = KS.f[ppm] - ffeq_ppm - Dsppm;
		const dreal Dhppz = KS.f[ppz] - ffeq_ppz - Dsppz;
		const dreal Dhppp = KS.f[ppp] - ffeq_ppp - Dsppp;

		const dreal beta = _beta;
		const dreal gamma = _gamma;

		// Forcing - BGK-like
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

		KS.f[mmm] -= beta*(no2*Dsmmm + gamma*Dhmmm) - (no1 - beta)*Smmm*ffeq_mmm;
		KS.f[mmz] -= beta*(no2*Dsmmz + gamma*Dhmmz) - (no1 - beta)*Smmz*ffeq_mmz;
		KS.f[mmp] -= beta*(no2*Dsmmp + gamma*Dhmmp) - (no1 - beta)*Smmp*ffeq_mmp;
		KS.f[mzm] -= beta*(no2*Dsmzm + gamma*Dhmzm) - (no1 - beta)*Smzm*ffeq_mzm;
		KS.f[mzz] -= beta*(no2*Dsmzz + gamma*Dhmzz) - (no1 - beta)*Smzz*ffeq_mzz;
		KS.f[mzp] -= beta*(no2*Dsmzp + gamma*Dhmzp) - (no1 - beta)*Smzp*ffeq_mzp;
		KS.f[mpm] -= beta*(no2*Dsmpm + gamma*Dhmpm) - (no1 - beta)*Smpm*ffeq_mpm;
		KS.f[mpz] -= beta*(no2*Dsmpz + gamma*Dhmpz) - (no1 - beta)*Smpz*ffeq_mpz;
		KS.f[mpp] -= beta*(no2*Dsmpp + gamma*Dhmpp) - (no1 - beta)*Smpp*ffeq_mpp;
		KS.f[zmm] -= beta*(no2*Dszmm + gamma*Dhzmm) - (no1 - beta)*Szmm*ffeq_zmm;
		KS.f[zmz] -= beta*(no2*Dszmz + gamma*Dhzmz) - (no1 - beta)*Szmz*ffeq_zmz;
		KS.f[zmp] -= beta*(no2*Dszmp + gamma*Dhzmp) - (no1 - beta)*Szmp*ffeq_zmp;
		KS.f[zzm] -= beta*(no2*Dszzm + gamma*Dhzzm) - (no1 - beta)*Szzm*ffeq_zzm;
		KS.f[zzz] -= beta*(no2*Dszzz + gamma*Dhzzz) - (no1 - beta)*Szzz*ffeq_zzz;
		KS.f[zzp] -= beta*(no2*Dszzp + gamma*Dhzzp) - (no1 - beta)*Szzp*ffeq_zzp;
		KS.f[zpm] -= beta*(no2*Dszpm + gamma*Dhzpm) - (no1 - beta)*Szpm*ffeq_zpm;
		KS.f[zpz] -= beta*(no2*Dszpz + gamma*Dhzpz) - (no1 - beta)*Szpz*ffeq_zpz;
		KS.f[zpp] -= beta*(no2*Dszpp + gamma*Dhzpp) - (no1 - beta)*Szpp*ffeq_zpp;
		KS.f[pmm] -= beta*(no2*Dspmm + gamma*Dhpmm) - (no1 - beta)*Spmm*ffeq_pmm;
		KS.f[pmz] -= beta*(no2*Dspmz + gamma*Dhpmz) - (no1 - beta)*Spmz*ffeq_pmz;
		KS.f[pmp] -= beta*(no2*Dspmp + gamma*Dhpmp) - (no1 - beta)*Spmp*ffeq_pmp;
		KS.f[pzm] -= beta*(no2*Dspzm + gamma*Dhpzm) - (no1 - beta)*Spzm*ffeq_pzm;
		KS.f[pzz] -= beta*(no2*Dspzz + gamma*Dhpzz) - (no1 - beta)*Spzz*ffeq_pzz;
		KS.f[pzp] -= beta*(no2*Dspzp + gamma*Dhpzp) - (no1 - beta)*Spzp*ffeq_pzp;
		KS.f[ppm] -= beta*(no2*Dsppm + gamma*Dhppm) - (no1 - beta)*Sppm*ffeq_ppm;
		KS.f[ppz] -= beta*(no2*Dsppz + gamma*Dhppz) - (no1 - beta)*Sppz*ffeq_ppz;
		KS.f[ppp] -= beta*(no2*Dsppp + gamma*Dhppp) - (no1 - beta)*Sppp*ffeq_ppp;
	}
};
