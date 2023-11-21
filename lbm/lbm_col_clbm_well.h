#include "lbm_common_well.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT_WELL<TRAITS>
>
struct LBM_CLBM_WELL : LBM_COMMON_WELL< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id =  "CLBM_WELL"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
#define K_mm0  n1o36
#define K_zm0  n1o9
#define K_pm0  n1o36
#define K_mz0  n1o9
#define K_zz0  n4o9
#define K_pz0  n1o9
#define K_mp0  n1o36
#define K_zp0  n1o9
#define K_pp0  n1o36

#define K_m00  n1o9
#define K_z00  n4o9
#define K_p00  n1o9
#define K_m01  0
#define K_z01  0
#define K_p01  0
#define K_m02  n1o27
#define K_z02  n4o27
#define K_p02  n1o27

#define K_000  1
#define K_010  0
#define K_020  n1o3
#define K_001  0
#define K_011  0
#define K_021  0
#define K_002  n1o3
#define K_012  0
#define K_022  n1o9

		// correction DEBUG:: set all "f" to well-conditioned ones
		// gen1.php BEGIN
		const dreal k_mm0 = (KS.f[mmp] + KS.f[mmm]) + KS.f[mmz];
		const dreal k_mz0 = (KS.f[mzp] + KS.f[mzm]) + KS.f[mzz];
		const dreal k_mp0 = (KS.f[mpp] + KS.f[mpm]) + KS.f[mpz];
		const dreal k_zm0 = (KS.f[zmp] + KS.f[zmm]) + KS.f[zmz];
		const dreal k_zz0 = (KS.f[zzp] + KS.f[zzm]) + KS.f[zzz];
		const dreal k_zp0 = (KS.f[zpp] + KS.f[zpm]) + KS.f[zpz];
		const dreal k_pm0 = (KS.f[pmp] + KS.f[pmm]) + KS.f[pmz];
		const dreal k_pz0 = (KS.f[pzp] + KS.f[pzm]) + KS.f[pzz];
		const dreal k_pp0 = (KS.f[ppp] + KS.f[ppm]) + KS.f[ppz];


		//Eq 7
		const dreal k_mm1 = (KS.f[mmp] - KS.f[mmm]) - KS.vz * ( k_mm0 + K_mm0);
		const dreal k_mz1 = (KS.f[mzp] - KS.f[mzm]) - KS.vz * ( k_mz0 + K_mz0);
		const dreal k_mp1 = (KS.f[mpp] - KS.f[mpm]) - KS.vz * ( k_mp0 + K_mp0);
		const dreal k_zm1 = (KS.f[zmp] - KS.f[zmm]) - KS.vz * ( k_zm0 + K_zm0);
		const dreal k_zz1 = (KS.f[zzp] - KS.f[zzm]) - KS.vz * ( k_zz0 + K_zz0);
		const dreal k_zp1 = (KS.f[zpp] - KS.f[zpm]) - KS.vz * ( k_zp0 + K_zp0);
		const dreal k_pm1 = (KS.f[pmp] - KS.f[pmm]) - KS.vz * ( k_pm0 + K_pm0);
		const dreal k_pz1 = (KS.f[pzp] - KS.f[pzm]) - KS.vz * ( k_pz0 + K_pz0);
		const dreal k_pp1 = (KS.f[ppp] - KS.f[ppm]) - KS.vz * ( k_pp0 + K_pp0);

		//Eq 8
		const dreal k_mm2 = (KS.f[mmp] + KS.f[mmm]) - no2*KS.vz*(KS.f[mmp] - KS.f[mmm]) + KS.vz*KS.vz*( k_mm0 + K_mm0);
		const dreal k_mz2 = (KS.f[mzp] + KS.f[mzm]) - no2*KS.vz*(KS.f[mzp] - KS.f[mzm]) + KS.vz*KS.vz*( k_mz0 + K_mz0);
		const dreal k_mp2 = (KS.f[mpp] + KS.f[mpm]) - no2*KS.vz*(KS.f[mpp] - KS.f[mpm]) + KS.vz*KS.vz*( k_mp0 + K_mp0);
		const dreal k_zm2 = (KS.f[zmp] + KS.f[zmm]) - no2*KS.vz*(KS.f[zmp] - KS.f[zmm]) + KS.vz*KS.vz*( k_zm0 + K_zm0);
		const dreal k_zz2 = (KS.f[zzp] + KS.f[zzm]) - no2*KS.vz*(KS.f[zzp] - KS.f[zzm]) + KS.vz*KS.vz*( k_zz0 + K_zz0);
		const dreal k_zp2 = (KS.f[zpp] + KS.f[zpm]) - no2*KS.vz*(KS.f[zpp] - KS.f[zpm]) + KS.vz*KS.vz*( k_zp0 + K_zp0);
		const dreal k_pm2 = (KS.f[pmp] + KS.f[pmm]) - no2*KS.vz*(KS.f[pmp] - KS.f[pmm]) + KS.vz*KS.vz*( k_pm0 + K_pm0);
		const dreal k_pz2 = (KS.f[pzp] + KS.f[pzm]) - no2*KS.vz*(KS.f[pzp] - KS.f[pzm]) + KS.vz*KS.vz*( k_pz0 + K_pz0);
		const dreal k_pp2 = (KS.f[ppp] + KS.f[ppm]) - no2*KS.vz*(KS.f[ppp] - KS.f[ppm]) + KS.vz*KS.vz*( k_pp0 + K_pp0);

		//Eq 9
		const dreal k_m00 = (k_mp0 + k_mm0) + k_mz0;
		const dreal k_z00 = (k_zp0 + k_zm0) + k_zz0;
		const dreal k_p00 = (k_pp0 + k_pm0) + k_pz0;
		const dreal k_m01 = (k_mp1 + k_mm1) + k_mz1;
		const dreal k_z01 = (k_zp1 + k_zm1) + k_zz1;
		const dreal k_p01 = (k_pp1 + k_pm1) + k_pz1;
		const dreal k_m02 = (k_mp2 + k_mm2) + k_mz2;
		const dreal k_z02 = (k_zp2 + k_zm2) + k_zz2;
		const dreal k_p02 = (k_pp2 + k_pm2) + k_pz2;

		//Eq 10
		const dreal k_m10 = (k_mp0 - k_mm0) - KS.vy * (  k_m00 + K_m00 );
		const dreal k_z10 = (k_zp0 - k_zm0) - KS.vy * (  k_z00 + K_z00 );
		const dreal k_p10 = (k_pp0 - k_pm0) - KS.vy * (  k_p00 + K_p00 );
		const dreal k_m11 = (k_mp1 - k_mm1) - KS.vy * (  k_m01 + K_m01 );
		const dreal k_z11 = (k_zp1 - k_zm1) - KS.vy * (  k_z01 + K_z01 );
		const dreal k_p11 = (k_pp1 - k_pm1) - KS.vy * (  k_p01 + K_p01 );
		const dreal k_m12 = (k_mp2 - k_mm2) - KS.vy * (  k_m02 + K_m02 );
		const dreal k_z12 = (k_zp2 - k_zm2) - KS.vy * (  k_z02 + K_z02 );
		const dreal k_p12 = (k_pp2 - k_pm2) - KS.vy * (  k_p02 + K_p02 );

		//Eq 11
		const dreal k_m20 = (k_mp0 + k_mm0) - no2*KS.vy* (k_mp0 - k_mm0) + KS.vy*KS.vy * (  k_m00 + K_m00 );
		const dreal k_z20 = (k_zp0 + k_zm0) - no2*KS.vy* (k_zp0 - k_zm0) + KS.vy*KS.vy * (  k_z00 + K_z00 );
		const dreal k_p20 = (k_pp0 + k_pm0) - no2*KS.vy* (k_pp0 - k_pm0) + KS.vy*KS.vy * (  k_p00 + K_p00 );
		const dreal k_m21 = (k_mp1 + k_mm1) - no2*KS.vy* (k_mp1 - k_mm1) + KS.vy*KS.vy * (  k_m01 + K_m01 );
		const dreal k_z21 = (k_zp1 + k_zm1) - no2*KS.vy* (k_zp1 - k_zm1) + KS.vy*KS.vy * (  k_z01 + K_z01 );
		const dreal k_p21 = (k_pp1 + k_pm1) - no2*KS.vy* (k_pp1 - k_pm1) + KS.vy*KS.vy * (  k_p01 + K_p01 );
		const dreal k_m22 = (k_mp2 + k_mm2) - no2*KS.vy* (k_mp2 - k_mm2) + KS.vy*KS.vy * (  k_m02 + K_m02 );
		const dreal k_z22 = (k_zp2 + k_zm2) - no2*KS.vy* (k_zp2 - k_zm2) + KS.vy*KS.vy * (  k_z02 + K_z02 );
		const dreal k_p22 = (k_pp2 + k_pm2) - no2*KS.vy* (k_pp2 - k_pm2) + KS.vy*KS.vy * (  k_p02 + K_p02 );

		//Eq 12
		dreal k_000 = (k_p00 + k_m00) + k_z00;
		dreal k_001 = (k_p01 + k_m01) + k_z01;
		dreal k_010 = (k_p10 + k_m10) + k_z10;
		dreal k_011 = (k_p11 + k_m11) + k_z11;
		dreal k_012 = (k_p12 + k_m12) + k_z12;
		dreal k_021 = (k_p21 + k_m21) + k_z21;
		dreal k_022 = (k_p22 + k_m22) + k_z22;
		dreal k_002 = (k_p02 + k_m02) + k_z02;
		dreal k_020 = (k_p20 + k_m20) + k_z20;

		//Eq 13
		dreal k_100 = (k_p00 - k_m00) - KS.vx * (  k_000 + K_000 );
		dreal k_101 = (k_p01 - k_m01) - KS.vx * (  k_001 + K_001 );
		dreal k_102 = (k_p02 - k_m02) - KS.vx * (  k_002 + K_002 );
		dreal k_110 = (k_p10 - k_m10) - KS.vx * (  k_010 + K_010 );
		dreal k_111 = (k_p11 - k_m11) - KS.vx * (  k_011 + K_011 );
		dreal k_112 = (k_p12 - k_m12) - KS.vx * (  k_012 + K_012 );
		dreal k_120 = (k_p20 - k_m20) - KS.vx * (  k_020 + K_020 );
		dreal k_121 = (k_p21 - k_m21) - KS.vx * (  k_021 + K_021 );
		dreal k_122 = (k_p22 - k_m22) - KS.vx * (  k_022 + K_022 );

		//Eq 14
		dreal k_201 = (k_p01 + k_m01) - no2*KS.vx* (k_p01 - k_m01) + KS.vx*KS.vx * (  k_001 + K_001 );
		dreal k_210 = (k_p10 + k_m10) - no2*KS.vx* (k_p10 - k_m10) + KS.vx*KS.vx * (  k_010 + K_010 );
		dreal k_211 = (k_p11 + k_m11) - no2*KS.vx* (k_p11 - k_m11) + KS.vx*KS.vx * (  k_011 + K_011 );
		dreal k_212 = (k_p12 + k_m12) - no2*KS.vx* (k_p12 - k_m12) + KS.vx*KS.vx * (  k_012 + K_012 );
		dreal k_221 = (k_p21 + k_m21) - no2*KS.vx* (k_p21 - k_m21) + KS.vx*KS.vx * (  k_021 + K_021 );
		dreal k_200 = (k_p00 + k_m00) - no2*KS.vx* (k_p00 - k_m00) + KS.vx*KS.vx * (  k_000 + K_000 );
		dreal k_202 = (k_p02 + k_m02) - no2*KS.vx* (k_p02 - k_m02) + KS.vx*KS.vx * (  k_002 + K_002 );
		dreal k_220 = (k_p20 + k_m20) - no2*KS.vx* (k_p20 - k_m20) + KS.vx*KS.vx * (  k_020 + K_020 );
		dreal k_222 = (k_p22 + k_m22) - no2*KS.vx* (k_p22 - k_m22) + KS.vx*KS.vx * (  k_022 + K_022 );

		// tak, tedka mam centralni momenty podle puvodniho clanku takove, ze se rovnaji tem no-well-cond
		// gen1.php END

		// relaxation definition
		dreal omega1 = no1/(no3*KS.lbmViscosity+n1o2);
		dreal omega2 = no1;//omega1;
		dreal omega3 = no1;
		dreal omega4 = no1;
		dreal omega5 = no1;
		dreal omega6 = no1;
		dreal omega7 = no1;
		dreal omega8 = no1;
		dreal omega9 = no1;
		dreal omega10 = no1;

		// derivatives of v: notation taken from Geier's paper 2015: Appendix D Eq D.1-3
		const dreal Dxu = - omega1/no2/KS.rho*(no2*k_200-k_020-k_002) - omega2/no2/KS.rho*(k_200+k_020+k_002 - k_000);
		const dreal Dyv = Dxu + n3o2*omega1/KS.rho *(k_200-k_020);
		const dreal Dzw = Dxu + n3o2*omega1/KS.rho *(k_200-k_002);

		// Eqs D.4-6
		const dreal EqD4RHS = (no1-omega1)*(k_200-k_020) - no3*KS.rho*(no1 - omega1*n1o2)*(KS.vx*KS.vx*Dxu - KS.vy*KS.vy*Dyv);
		const dreal EqD5RHS = (no1-omega1)*(k_200-k_002) - no3*KS.rho*(no1 - omega1*n1o2)*(KS.vx*KS.vx*Dxu - KS.vz*KS.vz*Dzw);
		const dreal EqD6RHS = k_000*omega2 + (no1-omega2)*(k_200+k_020+k_002) - no3*KS.rho*(no1-omega2/no2)*(KS.vx*KS.vx*Dxu + KS.vy*KS.vy*Dyv + KS.vz*KS.vz*Dzw);
		// see rovnice2.mw
		dreal ks_200 = n1o3*(EqD4RHS + EqD5RHS + EqD6RHS);
		dreal ks_020 = n1o3*(-no2*EqD4RHS + EqD5RHS + EqD6RHS);
		dreal ks_002 = n1o3*(EqD4RHS - no2*EqD5RHS + EqD6RHS);
		// Eqs D7-12: see lbm_cum.h
		const dreal ks_120 = (-k_102-k_120)*omega3*n1o2+(k_102-k_120)*omega4*n1o2+k_120;
		const dreal ks_102 = (-k_102-k_120)*omega3*n1o2+(-k_102+k_120)*omega4*n1o2+k_102;
		const dreal ks_210 = (-k_012-k_210)*omega3*n1o2+(k_012-k_210)*omega4*n1o2+k_210;
		const dreal ks_012 = (-k_012-k_210)*omega3*n1o2+(-k_012+k_210)*omega4*n1o2+k_012;
		const dreal ks_021 = (-k_021-k_201)*omega3*n1o2+(-k_021+k_201)*omega4*n1o2+k_021;
		const dreal ks_201 = (-k_021-k_201)*omega3*n1o2+(k_021-k_201)*omega4*n1o2+k_201;
		// Eq D13
		const dreal ks_111 = (no1-omega5)*k_111;
		// Eqs 43-45
		const dreal EqD14RHS = (no1-omega6)*(k_220-no2*k_202+k_022);
		const dreal EqD15RHS = (no1-omega6)*(k_220+k_202-no2*k_022);
		const dreal EqD16RHS = (no1-omega7)*(k_220+k_202+k_022) + n1o3*omega7*k_000; // FIXME
		// see rovnice2.mw
		dreal ks_220 = n1o3*(EqD14RHS+EqD15RHS+EqD16RHS);
		dreal ks_202 = n1o3*(-EqD14RHS+EqD16RHS);
		dreal ks_022 = n1o3*(-EqD15RHS+EqD16RHS);
		// Eq 46-48
		const dreal ks_211 = (no1-omega8)*k_211;
		const dreal ks_121 = (no1-omega8)*k_121;
		const dreal ks_112 = (no1-omega8)*k_112;
		// Eqs 49-52
		const dreal ks_221 = (no1-omega9)*k_221;
		const dreal ks_212 = (no1-omega9)*k_212;
		const dreal ks_122 = (no1-omega9)*k_122;
		dreal ks_222 = (no1-omega10)*k_222 + n1o27*omega10*k_000; // FIXME
//		const dreal ks_222 = (no1-omega10)*k_222;

		// backward central moment transformation
		// Geier 2017: forcing scheme ---> zde se musi dat (-1) pro forcing ... ale s nulovou silou a - zde to nefunguje ... 
		dreal ks_000 = k_000;
		const dreal ks_100 = k_100;
		const dreal ks_010 = k_010;
		const dreal ks_001 = k_001;
		
		// ad-hoc ?? FIXME
		const dreal ks_101 =  (no1 - omega1)*k_101;
		const dreal ks_011 =  (no1 - omega1)*k_011;
		const dreal ks_110 =  (no1 - omega1)*k_110;

		// gen3.php BEGIN
		//Eq57
		const dreal ks_z00 = ks_000*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_100 - ks_200 - K_000*KS.vx*KS.vx;
		const dreal ks_z01 = ks_001*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_101 - ks_201 - K_001*KS.vx*KS.vx;
		const dreal ks_z02 = ks_002*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_102 - ks_202 - K_002*KS.vx*KS.vx;
		const dreal ks_z10 = ks_010*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_110 - ks_210 - K_010*KS.vx*KS.vx;
		const dreal ks_z11 = ks_011*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_111 - ks_211 - K_011*KS.vx*KS.vx;
		const dreal ks_z12 = ks_012*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_112 - ks_212 - K_012*KS.vx*KS.vx;
		const dreal ks_z20 = ks_020*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_120 - ks_220 - K_020*KS.vx*KS.vx;
		const dreal ks_z21 = ks_021*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_121 - ks_221 - K_021*KS.vx*KS.vx;
		const dreal ks_z22 = ks_022*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_122 - ks_222 - K_022*KS.vx*KS.vx;

		//Eq58
		const dreal ks_m00 = ((ks_000 + K_000)*(KS.vx*KS.vx - KS.vx) + ks_100*(no2*KS.vx - no1) + ks_200)*n1o2;
		const dreal ks_m01 = ((ks_001 + K_001)*(KS.vx*KS.vx - KS.vx) + ks_101*(no2*KS.vx - no1) + ks_201)*n1o2;
		const dreal ks_m02 = ((ks_002 + K_002)*(KS.vx*KS.vx - KS.vx) + ks_102*(no2*KS.vx - no1) + ks_202)*n1o2;
		const dreal ks_m10 = ((ks_010 + K_010)*(KS.vx*KS.vx - KS.vx) + ks_110*(no2*KS.vx - no1) + ks_210)*n1o2;
		const dreal ks_m11 = ((ks_011 + K_011)*(KS.vx*KS.vx - KS.vx) + ks_111*(no2*KS.vx - no1) + ks_211)*n1o2;
		const dreal ks_m12 = ((ks_012 + K_012)*(KS.vx*KS.vx - KS.vx) + ks_112*(no2*KS.vx - no1) + ks_212)*n1o2;
		const dreal ks_m20 = ((ks_020 + K_020)*(KS.vx*KS.vx - KS.vx) + ks_120*(no2*KS.vx - no1) + ks_220)*n1o2;
		const dreal ks_m21 = ((ks_021 + K_021)*(KS.vx*KS.vx - KS.vx) + ks_121*(no2*KS.vx - no1) + ks_221)*n1o2;
		const dreal ks_m22 = ((ks_022 + K_022)*(KS.vx*KS.vx - KS.vx) + ks_122*(no2*KS.vx - no1) + ks_222)*n1o2;

		//Eq59
		const dreal ks_p00 = ((ks_000 + K_000)*(KS.vx*KS.vx + KS.vx) + ks_100*(no2*KS.vx + no1) + ks_200)*n1o2;
		const dreal ks_p01 = ((ks_001 + K_001)*(KS.vx*KS.vx + KS.vx) + ks_101*(no2*KS.vx + no1) + ks_201)*n1o2;
		const dreal ks_p02 = ((ks_002 + K_002)*(KS.vx*KS.vx + KS.vx) + ks_102*(no2*KS.vx + no1) + ks_202)*n1o2;
		const dreal ks_p10 = ((ks_010 + K_010)*(KS.vx*KS.vx + KS.vx) + ks_110*(no2*KS.vx + no1) + ks_210)*n1o2;
		const dreal ks_p11 = ((ks_011 + K_011)*(KS.vx*KS.vx + KS.vx) + ks_111*(no2*KS.vx + no1) + ks_211)*n1o2;
		const dreal ks_p12 = ((ks_012 + K_012)*(KS.vx*KS.vx + KS.vx) + ks_112*(no2*KS.vx + no1) + ks_212)*n1o2;
		const dreal ks_p20 = ((ks_020 + K_020)*(KS.vx*KS.vx + KS.vx) + ks_120*(no2*KS.vx + no1) + ks_220)*n1o2;
		const dreal ks_p21 = ((ks_021 + K_021)*(KS.vx*KS.vx + KS.vx) + ks_121*(no2*KS.vx + no1) + ks_221)*n1o2;
		const dreal ks_p22 = ((ks_022 + K_022)*(KS.vx*KS.vx + KS.vx) + ks_122*(no2*KS.vx + no1) + ks_222)*n1o2;

		//Eq60
		const dreal ks_mz0 = ks_m00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m10 - ks_m20 - K_m00*KS.vy*KS.vy;
		const dreal ks_mz1 = ks_m01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m11 - ks_m21 - K_m01*KS.vy*KS.vy;
		const dreal ks_mz2 = ks_m02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m12 - ks_m22 - K_m02*KS.vy*KS.vy;
		const dreal ks_zz0 = ks_z00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z10 - ks_z20 - K_z00*KS.vy*KS.vy;
		const dreal ks_zz1 = ks_z01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z11 - ks_z21 - K_z01*KS.vy*KS.vy;
		const dreal ks_zz2 = ks_z02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z12 - ks_z22 - K_z02*KS.vy*KS.vy;
		const dreal ks_pz0 = ks_p00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p10 - ks_p20 - K_p00*KS.vy*KS.vy;
		const dreal ks_pz1 = ks_p01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p11 - ks_p21 - K_p01*KS.vy*KS.vy;
		const dreal ks_pz2 = ks_p02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p12 - ks_p22 - K_p02*KS.vy*KS.vy;

		//Eq61
		const dreal ks_mm0 = ((ks_m00 + K_m00)*(KS.vy*KS.vy - KS.vy) + ks_m10*(no2*KS.vy - no1) + ks_m20)*n1o2;
		const dreal ks_mm1 = ((ks_m01 + K_m01)*(KS.vy*KS.vy - KS.vy) + ks_m11*(no2*KS.vy - no1) + ks_m21)*n1o2;
		const dreal ks_mm2 = ((ks_m02 + K_m02)*(KS.vy*KS.vy - KS.vy) + ks_m12*(no2*KS.vy - no1) + ks_m22)*n1o2;
		const dreal ks_zm0 = ((ks_z00 + K_z00)*(KS.vy*KS.vy - KS.vy) + ks_z10*(no2*KS.vy - no1) + ks_z20)*n1o2;
		const dreal ks_zm1 = ((ks_z01 + K_z01)*(KS.vy*KS.vy - KS.vy) + ks_z11*(no2*KS.vy - no1) + ks_z21)*n1o2;
		const dreal ks_zm2 = ((ks_z02 + K_z02)*(KS.vy*KS.vy - KS.vy) + ks_z12*(no2*KS.vy - no1) + ks_z22)*n1o2;
		const dreal ks_pm0 = ((ks_p00 + K_p00)*(KS.vy*KS.vy - KS.vy) + ks_p10*(no2*KS.vy - no1) + ks_p20)*n1o2;
		const dreal ks_pm1 = ((ks_p01 + K_p01)*(KS.vy*KS.vy - KS.vy) + ks_p11*(no2*KS.vy - no1) + ks_p21)*n1o2;
		const dreal ks_pm2 = ((ks_p02 + K_p02)*(KS.vy*KS.vy - KS.vy) + ks_p12*(no2*KS.vy - no1) + ks_p22)*n1o2;

		//Eq62
		const dreal ks_mp0 = ((ks_m00 + K_m00)*(KS.vy*KS.vy + KS.vy) + ks_m10*(no2*KS.vy + no1) + ks_m20)*n1o2;
		const dreal ks_mp1 = ((ks_m01 + K_m01)*(KS.vy*KS.vy + KS.vy) + ks_m11*(no2*KS.vy + no1) + ks_m21)*n1o2;
		const dreal ks_mp2 = ((ks_m02 + K_m02)*(KS.vy*KS.vy + KS.vy) + ks_m12*(no2*KS.vy + no1) + ks_m22)*n1o2;
		const dreal ks_zp0 = ((ks_z00 + K_z00)*(KS.vy*KS.vy + KS.vy) + ks_z10*(no2*KS.vy + no1) + ks_z20)*n1o2;
		const dreal ks_zp1 = ((ks_z01 + K_z01)*(KS.vy*KS.vy + KS.vy) + ks_z11*(no2*KS.vy + no1) + ks_z21)*n1o2;
		const dreal ks_zp2 = ((ks_z02 + K_z02)*(KS.vy*KS.vy + KS.vy) + ks_z12*(no2*KS.vy + no1) + ks_z22)*n1o2;
		const dreal ks_pp0 = ((ks_p00 + K_p00)*(KS.vy*KS.vy + KS.vy) + ks_p10*(no2*KS.vy + no1) + ks_p20)*n1o2;
		const dreal ks_pp1 = ((ks_p01 + K_p01)*(KS.vy*KS.vy + KS.vy) + ks_p11*(no2*KS.vy + no1) + ks_p21)*n1o2;
		const dreal ks_pp2 = ((ks_p02 + K_p02)*(KS.vy*KS.vy + KS.vy) + ks_p12*(no2*KS.vy + no1) + ks_p22)*n1o2;

		//Eq63
		KS.f[mmz] = ks_mm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mm1 - ks_mm2 - K_mm0*KS.vz*KS.vz;
		KS.f[mzz] = ks_mz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mz1 - ks_mz2 - K_mz0*KS.vz*KS.vz;
		KS.f[mpz] = ks_mp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mp1 - ks_mp2 - K_mp0*KS.vz*KS.vz;
		KS.f[zmz] = ks_zm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zm1 - ks_zm2 - K_zm0*KS.vz*KS.vz;
		KS.f[zzz] = ks_zz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zz1 - ks_zz2 - K_zz0*KS.vz*KS.vz;
		KS.f[zpz] = ks_zp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zp1 - ks_zp2 - K_zp0*KS.vz*KS.vz;
		KS.f[pmz] = ks_pm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pm1 - ks_pm2 - K_pm0*KS.vz*KS.vz;
		KS.f[pzz] = ks_pz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pz1 - ks_pz2 - K_pz0*KS.vz*KS.vz;
		KS.f[ppz] = ks_pp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pp1 - ks_pp2 - K_pp0*KS.vz*KS.vz;

		//Eq64
		KS.f[mmm] = ((ks_mm0 + K_mm0)*(KS.vz*KS.vz-KS.vz) + ks_mm1*(no2*KS.vz - no1) + ks_mm2)*n1o2;
		KS.f[mzm] = ((ks_mz0 + K_mz0)*(KS.vz*KS.vz-KS.vz) + ks_mz1*(no2*KS.vz - no1) + ks_mz2)*n1o2;
		KS.f[mpm] = ((ks_mp0 + K_mp0)*(KS.vz*KS.vz-KS.vz) + ks_mp1*(no2*KS.vz - no1) + ks_mp2)*n1o2;
		KS.f[zmm] = ((ks_zm0 + K_zm0)*(KS.vz*KS.vz-KS.vz) + ks_zm1*(no2*KS.vz - no1) + ks_zm2)*n1o2;
		KS.f[zzm] = ((ks_zz0 + K_zz0)*(KS.vz*KS.vz-KS.vz) + ks_zz1*(no2*KS.vz - no1) + ks_zz2)*n1o2;
		KS.f[zpm] = ((ks_zp0 + K_zp0)*(KS.vz*KS.vz-KS.vz) + ks_zp1*(no2*KS.vz - no1) + ks_zp2)*n1o2;
		KS.f[pmm] = ((ks_pm0 + K_pm0)*(KS.vz*KS.vz-KS.vz) + ks_pm1*(no2*KS.vz - no1) + ks_pm2)*n1o2;
		KS.f[pzm] = ((ks_pz0 + K_pz0)*(KS.vz*KS.vz-KS.vz) + ks_pz1*(no2*KS.vz - no1) + ks_pz2)*n1o2;
		KS.f[ppm] = ((ks_pp0 + K_pp0)*(KS.vz*KS.vz-KS.vz) + ks_pp1*(no2*KS.vz - no1) + ks_pp2)*n1o2;

		//Eq65
		KS.f[mmp] = ((ks_mm0 + K_mm0)*(KS.vz*KS.vz+KS.vz) + ks_mm1*(no2*KS.vz + no1) + ks_mm2)*n1o2;
		KS.f[mzp] = ((ks_mz0 + K_mz0)*(KS.vz*KS.vz+KS.vz) + ks_mz1*(no2*KS.vz + no1) + ks_mz2)*n1o2;
		KS.f[mpp] = ((ks_mp0 + K_mp0)*(KS.vz*KS.vz+KS.vz) + ks_mp1*(no2*KS.vz + no1) + ks_mp2)*n1o2;
		KS.f[zmp] = ((ks_zm0 + K_zm0)*(KS.vz*KS.vz+KS.vz) + ks_zm1*(no2*KS.vz + no1) + ks_zm2)*n1o2;
		KS.f[zzp] = ((ks_zz0 + K_zz0)*(KS.vz*KS.vz+KS.vz) + ks_zz1*(no2*KS.vz + no1) + ks_zz2)*n1o2;
		KS.f[zpp] = ((ks_zp0 + K_zp0)*(KS.vz*KS.vz+KS.vz) + ks_zp1*(no2*KS.vz + no1) + ks_zp2)*n1o2;
		KS.f[pmp] = ((ks_pm0 + K_pm0)*(KS.vz*KS.vz+KS.vz) + ks_pm1*(no2*KS.vz + no1) + ks_pm2)*n1o2;
		KS.f[pzp] = ((ks_pz0 + K_pz0)*(KS.vz*KS.vz+KS.vz) + ks_pz1*(no2*KS.vz + no1) + ks_pz2)*n1o2;
		KS.f[ppp] = ((ks_pp0 + K_pp0)*(KS.vz*KS.vz+KS.vz) + ks_pp1*(no2*KS.vz + no1) + ks_pp2)*n1o2;

		#undef K_mm0
		#undef K_zm0
		#undef K_pm0
		#undef K_mz0
		#undef K_zz0
		#undef K_pz0
		#undef K_mp0
		#undef K_zp0
		#undef K_pp0

		#undef K_m00 
		#undef K_z00 
		#undef K_p00 
		#undef K_m01 
		#undef K_z01 
		#undef K_p01 
		#undef K_m02 
		#undef K_z02 
		#undef K_p02 

		#undef K_000 
		#undef K_010 
		#undef K_020 
		#undef K_001 
		#undef K_011 
		#undef K_021 
		#undef K_002 
		#undef K_012 
		#undef K_022 
	}
};

