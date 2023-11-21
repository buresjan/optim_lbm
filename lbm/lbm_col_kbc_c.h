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
// C1-C4 models with the central moment representation of DFs

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
//central moments - models C1, C2, C3, C4
// special moments with physical meaning & their equilibria
#define cT    (rT - KS.rho*(KS.vx*KS.vx + KS.vy*KS.vy + KS.vz*KS.vz))
#define cNxz  (rNxz + KS.rho*(KS.vz*KS.vz - KS.vx*KS.vx))
#define cNyz  (rNyz + KS.rho*(KS.vz*KS.vz - KS.vy*KS.vy))
#define cPxy  (rPxy - KS.rho*KS.vx*KS.vy)
#define cPxz  (rPxz - KS.rho*KS.vx*KS.vz)
#define cPyz  (rPyz - KS.rho*KS.vy*KS.vz)
#define cQxxy (rQxxy - n1o3*(no6*KS.vx*cPxy + KS.vy*(no3*KS.vx*KS.vx + no2*cNxz - cNyz + cT)))
#define cQxxz (rQxxz - n1o3*(no6*KS.vx*cPxz + KS.vz*(no3*KS.vx*KS.vx + no2*cNxz - cNyz + cT)))
#define cQxyy (rQxyy - n1o3*(no6*KS.vy*cPxy + KS.vx*(no3*KS.vy*KS.vy + no2*cNyz - cNxz + cT)))
#define cQyyz (rQyyz - n1o3*(no6*KS.vy*cPyz + KS.vz*(no3*KS.vy*KS.vy + no2*cNyz - cNxz + cT)))
#define cQxzz (rQxzz - n1o3*(no6*KS.vz*cPxz + KS.vx*(no3*KS.vz*KS.vz - cNyz - cNxz + cT)))
#define cQyzz (rQyzz - n1o3*(no6*KS.vz*cPyz + KS.vy*(no3*KS.vz*KS.vz - cNyz - cNxz + cT)))
#define cQxyz (rQxyz - KS.vx*cPyz - KS.vy*cPxz - KS.vz*cPxy - KS.vx*KS.vy*KS.vz)
//  their equilibria
#define ecT   (KS.rho*no3*_csqr)
#define ecNxz no0
#define ecNyz no0
#define ecPxy no0
#define ecPxz no0
#define ecPyz no0
#define ecQxxy no0
#define ecQxxz no0
#define ecQxyy no0
#define ecQyyz no0
#define ecQxzz no0
#define ecQyzz no0
#define ecQxyz no0
// tensors included in the shear part FIXME -- this is not so simple, there must be prefactors with vx,vy,vz see Appendix in the KBC paper
//  deviatoric stress tensor D
#define cDzzz no0
#define cDpzz ((no2*cNxz - cNyz)/no6)
#define cDmzz ((no2*cNxz - cNyz)/no6)
#define cDzpz ((-cNxz + no2*cNyz)/no6)
#define cDzmz ((-cNxz + no2*cNyz)/no6)
#define cDzzp ((-cNxz - cNyz)/no6)
#define cDzzm ((-cNxz - cNyz)/no6)
#define cDmmz (cPxy/no4)
#define cDmpz (-cPxy/no4)
#define cDpmz (-cPxy/no4)
#define cDppz (cPxy/no4)
#define cDmzm (cPxz/no4)
#define cDmzp (-cPxz/no4)
#define cDpzm (-cPxz/no4)
#define cDpzp (cPxz/no4)
#define cDzmm (cPyz/no4)
#define cDzmp (-cPyz/no4)
#define cDzpm (-cPyz/no4)
#define cDzpp (cPyz/no4)
#define cDmmm no0
#define cDmmp no0
#define cDmpm no0
#define cDmpp no0
#define cDpmm no0
#define cDpmp no0
#define cDppm no0
#define cDppp no0

#define ecDzzz no0
#define ecDpzz no0
#define ecDmzz no0
#define ecDzpz no0
#define ecDzmz no0
#define ecDzzp no0
#define ecDzzm no0
#define ecDmmz no0
#define ecDmpz no0
#define ecDpmz no0
#define ecDppz no0
#define ecDmzm no0
#define ecDmzp no0
#define ecDpzm no0
#define ecDpzp no0
#define ecDzmm no0
#define ecDzmp no0
#define ecDzpm no0
#define ecDzpp no0
#define ecDmmm no0
#define ecDmmp no0
#define ecDmpm no0
#define ecDmpp no0
#define ecDpmm no0
#define ecDpmp no0
#define ecDppm no0
#define ecDppp no0
//  trace of stress tensor T
#define cTzzz -cT
#define cTpzz (cT/no6)
#define cTmzz (cT/no6)
#define cTzpz (cT/no6)
#define cTzmz (cT/no6)
#define cTzzp (cT/no6)
#define cTzzm (cT/no6)
#define cTmmz no0
#define cTmpz no0
#define cTpmz no0
#define cTppz no0
#define cTmzm no0
#define cTmzp no0
#define cTpzm no0
#define cTpzp no0
#define cTzmm no0
#define cTzmp no0
#define cTzpm no0
#define cTzpp no0
#define cTmmm no0
#define cTmmp no0
#define cTmpm no0
#define cTmpp no0
#define cTpmm no0
#define cTpmp no0
#define cTppm no0
#define cTppp no0

#define ecTzzz -ecT
#define ecTpzz (ecT/no6)
#define ecTmzz (ecT/no6)
#define ecTzpz (ecT/no6)
#define ecTzmz (ecT/no6)
#define ecTzzp (ecT/no6)
#define ecTzzm (ecT/no6)
#define ecTmmz no0
#define ecTmpz no0
#define ecTpmz no0
#define ecTppz no0
#define ecTmzm no0
#define ecTmzp no0
#define ecTpzm no0
#define ecTpzp no0
#define ecTzmm no0
#define ecTzmp no0
#define ecTzpm no0
#define ecTzpp no0
#define ecTmmm no0
#define ecTmmp no0
#define ecTmpm no0
#define ecTmpp no0
#define ecTpmm no0
#define ecTpmp no0
#define ecTppm no0
#define ecTppp no0
//  heat flux tensor Q
#define cQzzz no0
#define cQpzz (-(cQxyy + cQxzz)/no2)
#define cQmzz (+(cQxyy + cQxzz)/no2)
#define cQzpz (-(cQxxy + cQyzz)/no2)
#define cQzmz (+(cQxxy + cQyzz)/no2)
#define cQzzp (-(cQxxz + cQyyz)/no2)
#define cQzzm (+(cQxxz + cQyyz)/no2)
#define cQmmz ((-cQxyy - cQxxy)/no4)
#define cQmpz ((-cQxyy + cQxxy)/no4)
#define cQpmz ((+cQxyy - cQxxy)/no4)
#define cQppz ((+cQxyy + cQxxy)/no4)
#define cQmzm ((-cQxzz - cQxxz)/no4)
#define cQmzp ((-cQxzz + cQxxz)/no4)
#define cQpzm ((+cQxzz - cQxxz)/no4)
#define cQpzp ((+cQxzz + cQxxz)/no4)
#define cQzmm ((-cQyzz - cQyyz)/no4)
#define cQzmp ((-cQyzz + cQyyz)/no4)
#define cQzpm ((+cQyzz - cQyyz)/no4)
#define cQzpp ((+cQyzz + cQyyz)/no4)
#define cQmmm (-cQxyz/no8)
#define cQmmp (+cQxyz/no8)
#define cQmpm (+cQxyz/no8)
#define cQmpp (-cQxyz/no8)
#define cQpmm (+cQxyz/no8)
#define cQpmp (-cQxyz/no8)
#define cQppm (-cQxyz/no8)
#define cQppp (+cQxyz/no8)

#define ecQzzz no0
#define ecQpzz no0
#define ecQmzz no0
#define ecQzpz no0
#define ecQzmz no0
#define ecQzzp no0
#define ecQzzm no0
#define ecQmmz no0
#define ecQmpz no0
#define ecQpmz no0
#define ecQppz no0
#define ecQmzm no0
#define ecQmzp no0
#define ecQpzm no0
#define ecQpzp no0
#define ecQzmm no0
#define ecQzmp no0
#define ecQzpm no0
#define ecQzpp no0
#define ecQmmm no0
#define ecQmmp no0
#define ecQmpm no0
#define ecQmpp no0
#define ecQpmm no0
#define ecQpmp no0
#define ecQppm no0
#define ecQppp no0

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
struct LBM_KBC_C1 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_C1";
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

		const dreal ifeq_mmm = no1/ffeq_mmm;
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

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D~
		const dreal	Dsmmm = cDmmm - ecDmmm;
		const dreal Dsmmz = cDmmz - ecDmmz;
		const dreal Dsmmp = cDmmp - ecDmmp;
		const dreal Dsmzm = cDmzm - ecDmzm;
		const dreal Dsmzz = cDmzz - ecDmzz;
		const dreal Dsmzp = cDmzp - ecDmzp;
		const dreal Dsmpm = cDmpm - ecDmpm;
		const dreal Dsmpz = cDmpz - ecDmpz;
		const dreal Dsmpp = cDmpp - ecDmpp;
		const dreal Dszmm = cDzmm - ecDzmm;
		const dreal Dszmz = cDzmz - ecDzmz;
		const dreal Dszmp = cDzmp - ecDzmp;
		const dreal Dszzm = cDzzm - ecDzzm;
		const dreal Dszzz = cDzzz - ecDzzz;
		const dreal Dszzp = cDzzp - ecDzzp;
		const dreal Dszpm = cDzpm - ecDzpm;
		const dreal Dszpz = cDzpz - ecDzpz;
		const dreal Dszpp = cDzpp - ecDzpp;
		const dreal Dspmm = cDpmm - ecDpmm;
		const dreal Dspmz = cDpmz - ecDpmz;
		const dreal Dspmp = cDpmp - ecDpmp;
		const dreal Dspzm = cDpzm - ecDpzm;
		const dreal Dspzz = cDpzz - ecDpzz;
		const dreal Dspzp = cDpzp - ecDpzp;
		const dreal Dsppm = cDppm - ecDppm;
		const dreal Dsppz = cDppz - ecDppz;
		const dreal Dsppp = cDppp - ecDppp;

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
struct LBM_KBC_C2 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_C2";
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

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D~ + T~
		const dreal	Dsmmm = cDmmm - ecDmmm + cTmmm - ecTmmm;
		const dreal Dsmmz = cDmmz - ecDmmz + cTmmz - ecTmmz;
		const dreal Dsmmp = cDmmp - ecDmmp + cTmmp - ecTmmp;
		const dreal Dsmzm = cDmzm - ecDmzm + cTmzm - ecTmzm;
		const dreal Dsmzz = cDmzz - ecDmzz + cTmzz - ecTmzz;
		const dreal Dsmzp = cDmzp - ecDmzp + cTmzp - ecTmzp;
		const dreal Dsmpm = cDmpm - ecDmpm + cTmpm - ecTmpm;
		const dreal Dsmpz = cDmpz - ecDmpz + cTmpz - ecTmpz;
		const dreal Dsmpp = cDmpp - ecDmpp + cTmpp - ecTmpp;
		const dreal Dszmm = cDzmm - ecDzmm + cTzmm - ecTzmm;
		const dreal Dszmz = cDzmz - ecDzmz + cTzmz - ecTzmz;
		const dreal Dszmp = cDzmp - ecDzmp + cTzmp - ecTzmp;
		const dreal Dszzm = cDzzm - ecDzzm + cTzzm - ecTzzm;
		const dreal Dszzz = cDzzz - ecDzzz + cTzzz - ecTzzz;
		const dreal Dszzp = cDzzp - ecDzzp + cTzzp - ecTzzp;
		const dreal Dszpm = cDzpm - ecDzpm + cTzpm - ecTzpm;
		const dreal Dszpz = cDzpz - ecDzpz + cTzpz - ecTzpz;
		const dreal Dszpp = cDzpp - ecDzpp + cTzpp - ecTzpp;
		const dreal Dspmm = cDpmm - ecDpmm + cTpmm - ecTpmm;
		const dreal Dspmz = cDpmz - ecDpmz + cTpmz - ecTpmz;
		const dreal Dspmp = cDpmp - ecDpmp + cTpmp - ecTpmp;
		const dreal Dspzm = cDpzm - ecDpzm + cTpzm - ecTpzm;
		const dreal Dspzz = cDpzz - ecDpzz + cTpzz - ecTpzz;
		const dreal Dspzp = cDpzp - ecDpzp + cTpzp - ecTpzp;
		const dreal Dsppm = cDppm - ecDppm + cTppm - ecTppm;
		const dreal Dsppz = cDppz - ecDppz + cTppz - ecTppz;
		const dreal Dsppp = cDppp - ecDppp + cTppp - ecTppp;

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
struct LBM_KBC_C3 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_C3";
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

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D~ + Q~
		const dreal	Dsmmm = cDmmm - ecDmmm + cQmmm - ecQmmm;
		const dreal Dsmmz = cDmmz - ecDmmz + cQmmz - ecQmmz;
		const dreal Dsmmp = cDmmp - ecDmmp + cQmmp - ecQmmp;
		const dreal Dsmzm = cDmzm - ecDmzm + cQmzm - ecQmzm;
		const dreal Dsmzz = cDmzz - ecDmzz + cQmzz - ecQmzz;
		const dreal Dsmzp = cDmzp - ecDmzp + cQmzp - ecQmzp;
		const dreal Dsmpm = cDmpm - ecDmpm + cQmpm - ecQmpm;
		const dreal Dsmpz = cDmpz - ecDmpz + cQmpz - ecQmpz;
		const dreal Dsmpp = cDmpp - ecDmpp + cQmpp - ecQmpp;
		const dreal Dszmm = cDzmm - ecDzmm + cQzmm - ecQzmm;
		const dreal Dszmz = cDzmz - ecDzmz + cQzmz - ecQzmz;
		const dreal Dszmp = cDzmp - ecDzmp + cQzmp - ecQzmp;
		const dreal Dszzm = cDzzm - ecDzzm + cQzzm - ecQzzm;
		const dreal Dszzz = cDzzz - ecDzzz + cQzzz - ecQzzz;
		const dreal Dszzp = cDzzp - ecDzzp + cQzzp - ecQzzp;
		const dreal Dszpm = cDzpm - ecDzpm + cQzpm - ecQzpm;
		const dreal Dszpz = cDzpz - ecDzpz + cQzpz - ecQzpz;
		const dreal Dszpp = cDzpp - ecDzpp + cQzpp - ecQzpp;
		const dreal Dspmm = cDpmm - ecDpmm + cQpmm - ecQpmm;
		const dreal Dspmz = cDpmz - ecDpmz + cQpmz - ecQpmz;
		const dreal Dspmp = cDpmp - ecDpmp + cQpmp - ecQpmp;
		const dreal Dspzm = cDpzm - ecDpzm + cQpzm - ecQpzm;
		const dreal Dspzz = cDpzz - ecDpzz + cQpzz - ecQpzz;
		const dreal Dspzp = cDpzp - ecDpzp + cQpzp - ecQpzp;
		const dreal Dsppm = cDppm - ecDppm + cQppm - ecQppm;
		const dreal Dsppz = cDppz - ecDppz + cQppz - ecQppz;
		const dreal Dsppp = cDppp - ecDppp + cQppp - ecQppp;

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
struct LBM_KBC_C4 : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "KBC_C4";
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

		const dreal ifeq_mmm = no1/ffeq_mmm;
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

		// delta of shear part: delta s_i = s_i - s_i^eq; s = D~ + T~ + Q~
		const dreal	Dsmmm = cDmmm - ecDmmm + cTmmm - ecTmmm + cQmmm - ecQmmm;
		const dreal Dsmmz = cDmmz - ecDmmz + cTmmz - ecTmmz + cQmmz - ecQmmz;
		const dreal Dsmmp = cDmmp - ecDmmp + cTmmp - ecTmmp + cQmmp - ecQmmp;
		const dreal Dsmzm = cDmzm - ecDmzm + cTmzm - ecTmzm + cQmzm - ecQmzm;
		const dreal Dsmzz = cDmzz - ecDmzz + cTmzz - ecTmzz + cQmzz - ecQmzz;
		const dreal Dsmzp = cDmzp - ecDmzp + cTmzp - ecTmzp + cQmzp - ecQmzp;
		const dreal Dsmpm = cDmpm - ecDmpm + cTmpm - ecTmpm + cQmpm - ecQmpm;
		const dreal Dsmpz = cDmpz - ecDmpz + cTmpz - ecTmpz + cQmpz - ecQmpz;
		const dreal Dsmpp = cDmpp - ecDmpp + cTmpp - ecTmpp + cQmpp - ecQmpp;
		const dreal Dszmm = cDzmm - ecDzmm + cTzmm - ecTzmm + cQzmm - ecQzmm;
		const dreal Dszmz = cDzmz - ecDzmz + cTzmz - ecTzmz + cQzmz - ecQzmz;
		const dreal Dszmp = cDzmp - ecDzmp + cTzmp - ecTzmp + cQzmp - ecQzmp;
		const dreal Dszzm = cDzzm - ecDzzm + cTzzm - ecTzzm + cQzzm - ecQzzm;
		const dreal Dszzz = cDzzz - ecDzzz + cTzzz - ecTzzz + cQzzz - ecQzzz;
		const dreal Dszzp = cDzzp - ecDzzp + cTzzp - ecTzzp + cQzzp - ecQzzp;
		const dreal Dszpm = cDzpm - ecDzpm + cTzpm - ecTzpm + cQzpm - ecQzpm;
		const dreal Dszpz = cDzpz - ecDzpz + cTzpz - ecTzpz + cQzpz - ecQzpz;
		const dreal Dszpp = cDzpp - ecDzpp + cTzpp - ecTzpp + cQzpp - ecQzpp;
		const dreal Dspmm = cDpmm - ecDpmm + cTpmm - ecTpmm + cQpmm - ecQpmm;
		const dreal Dspmz = cDpmz - ecDpmz + cTpmz - ecTpmz + cQpmz - ecQpmz;
		const dreal Dspmp = cDpmp - ecDpmp + cTpmp - ecTpmp + cQpmp - ecQpmp;
		const dreal Dspzm = cDpzm - ecDpzm + cTpzm - ecTpzm + cQpzm - ecQpzm;
		const dreal Dspzz = cDpzz - ecDpzz + cTpzz - ecTpzz + cQpzz - ecQpzz;
		const dreal Dspzp = cDpzp - ecDpzp + cTpzp - ecTpzp + cQpzp - ecQpzp;
		const dreal Dsppm = cDppm - ecDppm + cTppm - ecTppm + cQppm - ecQppm;
		const dreal Dsppz = cDppz - ecDppz + cTppz - ecTppz + cQppz - ecQppz;
		const dreal Dsppp = cDppp - ecDppp + cTppp - ecTppp + cQppp - ecQppp;

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

