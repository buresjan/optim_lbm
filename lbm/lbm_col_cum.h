#include "lbm_common.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_CUM : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "CUM"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		// for compatibility of notation
		#define C_011 k_011
		#define C_010 k_010
		#define C_100 k_100
		#define C_101 k_101
		#define C_110 k_110
		#define C_111 k_111
		#define C_002 k_002
		#define C_020 k_020
		#define C_200 k_200
		#define C_012 k_012
		#define C_021 k_021
		#define C_102 k_102
		#define C_201 k_201
		#define C_120 k_120
		#define C_210 k_210
		// backward translation
		#define ks_001 Cs_001
		#define ks_011 Cs_011
		#define ks_010 Cs_010
		#define ks_100 Cs_100
		#define ks_101 Cs_101
		#define ks_110 Cs_110
		#define ks_111 Cs_111
		#define ks_002 Cs_002
		#define ks_020 Cs_020
		#define ks_200 Cs_200
		#define ks_012 Cs_012
		#define ks_021 Cs_021
		#define ks_102 Cs_102
		#define ks_201 Cs_201
		#define ks_120 Cs_120
		#define ks_210 Cs_210
		// backward translation
		// gen1nowell.php BEGIN
		//Eq 6
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
		const dreal k_mm1 = (KS.f[mmp] - KS.f[mmm]) - KS.vz * k_mm0;
		const dreal k_mz1 = (KS.f[mzp] - KS.f[mzm]) - KS.vz * k_mz0;
		const dreal k_mp1 = (KS.f[mpp] - KS.f[mpm]) - KS.vz * k_mp0;
		const dreal k_zm1 = (KS.f[zmp] - KS.f[zmm]) - KS.vz * k_zm0;
		const dreal k_zz1 = (KS.f[zzp] - KS.f[zzm]) - KS.vz * k_zz0;
		const dreal k_zp1 = (KS.f[zpp] - KS.f[zpm]) - KS.vz * k_zp0;
		const dreal k_pm1 = (KS.f[pmp] - KS.f[pmm]) - KS.vz * k_pm0;
		const dreal k_pz1 = (KS.f[pzp] - KS.f[pzm]) - KS.vz * k_pz0;
		const dreal k_pp1 = (KS.f[ppp] - KS.f[ppm]) - KS.vz * k_pp0;

		//Eq 8
		const dreal k_mm2 = (KS.f[mmp] + KS.f[mmm]) - no2*KS.vz*(KS.f[mmp] - KS.f[mmm]) + KS.vz*KS.vz* k_mm0;
		const dreal k_mz2 = (KS.f[mzp] + KS.f[mzm]) - no2*KS.vz*(KS.f[mzp] - KS.f[mzm]) + KS.vz*KS.vz* k_mz0;
		const dreal k_mp2 = (KS.f[mpp] + KS.f[mpm]) - no2*KS.vz*(KS.f[mpp] - KS.f[mpm]) + KS.vz*KS.vz* k_mp0;
		const dreal k_zm2 = (KS.f[zmp] + KS.f[zmm]) - no2*KS.vz*(KS.f[zmp] - KS.f[zmm]) + KS.vz*KS.vz* k_zm0;
		const dreal k_zz2 = (KS.f[zzp] + KS.f[zzm]) - no2*KS.vz*(KS.f[zzp] - KS.f[zzm]) + KS.vz*KS.vz* k_zz0;
		const dreal k_zp2 = (KS.f[zpp] + KS.f[zpm]) - no2*KS.vz*(KS.f[zpp] - KS.f[zpm]) + KS.vz*KS.vz* k_zp0;
		const dreal k_pm2 = (KS.f[pmp] + KS.f[pmm]) - no2*KS.vz*(KS.f[pmp] - KS.f[pmm]) + KS.vz*KS.vz* k_pm0;
		const dreal k_pz2 = (KS.f[pzp] + KS.f[pzm]) - no2*KS.vz*(KS.f[pzp] - KS.f[pzm]) + KS.vz*KS.vz* k_pz0;
		const dreal k_pp2 = (KS.f[ppp] + KS.f[ppm]) - no2*KS.vz*(KS.f[ppp] - KS.f[ppm]) + KS.vz*KS.vz* k_pp0;

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
		const dreal k_m10 = (k_mp0 - k_mm0) - KS.vy * k_m00;
		const dreal k_z10 = (k_zp0 - k_zm0) - KS.vy * k_z00;
		const dreal k_p10 = (k_pp0 - k_pm0) - KS.vy * k_p00;
		const dreal k_m11 = (k_mp1 - k_mm1) - KS.vy * k_m01;
		const dreal k_z11 = (k_zp1 - k_zm1) - KS.vy * k_z01;
		const dreal k_p11 = (k_pp1 - k_pm1) - KS.vy * k_p01;
		const dreal k_m12 = (k_mp2 - k_mm2) - KS.vy * k_m02;
		const dreal k_z12 = (k_zp2 - k_zm2) - KS.vy * k_z02;
		const dreal k_p12 = (k_pp2 - k_pm2) - KS.vy * k_p02;

		//Eq 11
		const dreal k_m20 = (k_mp0 + k_mm0) - no2*KS.vy* (k_mp0 - k_mm0) + KS.vy*KS.vy * k_m00;
		const dreal k_z20 = (k_zp0 + k_zm0) - no2*KS.vy* (k_zp0 - k_zm0) + KS.vy*KS.vy * k_z00;
		const dreal k_p20 = (k_pp0 + k_pm0) - no2*KS.vy* (k_pp0 - k_pm0) + KS.vy*KS.vy * k_p00;
		const dreal k_m21 = (k_mp1 + k_mm1) - no2*KS.vy* (k_mp1 - k_mm1) + KS.vy*KS.vy * k_m01;
		const dreal k_z21 = (k_zp1 + k_zm1) - no2*KS.vy* (k_zp1 - k_zm1) + KS.vy*KS.vy * k_z01;
		const dreal k_p21 = (k_pp1 + k_pm1) - no2*KS.vy* (k_pp1 - k_pm1) + KS.vy*KS.vy * k_p01;
		const dreal k_m22 = (k_mp2 + k_mm2) - no2*KS.vy* (k_mp2 - k_mm2) + KS.vy*KS.vy * k_m02;
		const dreal k_z22 = (k_zp2 + k_zm2) - no2*KS.vy* (k_zp2 - k_zm2) + KS.vy*KS.vy * k_z02;
		const dreal k_p22 = (k_pp2 + k_pm2) - no2*KS.vy* (k_pp2 - k_pm2) + KS.vy*KS.vy * k_p02;

		//Eq 12
		const dreal k_000 = (k_p00 + k_m00) + k_z00;
		const dreal k_001 = (k_p01 + k_m01) + k_z01;
		const dreal k_002 = (k_p02 + k_m02) + k_z02;
		const dreal k_010 = (k_p10 + k_m10) + k_z10;
		const dreal k_011 = (k_p11 + k_m11) + k_z11;
		const dreal k_012 = (k_p12 + k_m12) + k_z12;
		const dreal k_020 = (k_p20 + k_m20) + k_z20;
		const dreal k_021 = (k_p21 + k_m21) + k_z21;
		const dreal k_022 = (k_p22 + k_m22) + k_z22;

		//Eq 13
		const dreal k_100 = (k_p00 - k_m00) - KS.vx * k_000 ;
		const dreal k_101 = (k_p01 - k_m01) - KS.vx * k_001 ;
		const dreal k_102 = (k_p02 - k_m02) - KS.vx * k_002 ;
		const dreal k_110 = (k_p10 - k_m10) - KS.vx * k_010 ;
		const dreal k_111 = (k_p11 - k_m11) - KS.vx * k_011 ;
		const dreal k_112 = (k_p12 - k_m12) - KS.vx * k_012 ;
		const dreal k_120 = (k_p20 - k_m20) - KS.vx * k_020 ;
		const dreal k_121 = (k_p21 - k_m21) - KS.vx * k_021 ;
		const dreal k_122 = (k_p22 - k_m22) - KS.vx * k_022 ;

		//Eq 14
		const dreal k_200 = (k_p00 + k_m00) - no2*KS.vx* (k_p00 - k_m00) + KS.vx*KS.vx * k_000 ;
		const dreal k_201 = (k_p01 + k_m01) - no2*KS.vx* (k_p01 - k_m01) + KS.vx*KS.vx * k_001 ;
		const dreal k_202 = (k_p02 + k_m02) - no2*KS.vx* (k_p02 - k_m02) + KS.vx*KS.vx * k_002 ;
		const dreal k_210 = (k_p10 + k_m10) - no2*KS.vx* (k_p10 - k_m10) + KS.vx*KS.vx * k_010 ;
		const dreal k_211 = (k_p11 + k_m11) - no2*KS.vx* (k_p11 - k_m11) + KS.vx*KS.vx * k_011 ;
		const dreal k_212 = (k_p12 + k_m12) - no2*KS.vx* (k_p12 - k_m12) + KS.vx*KS.vx * k_012 ;
		const dreal k_220 = (k_p20 + k_m20) - no2*KS.vx* (k_p20 - k_m20) + KS.vx*KS.vx * k_020 ;
		const dreal k_221 = (k_p21 + k_m21) - no2*KS.vx* (k_p21 - k_m21) + KS.vx*KS.vx * k_021 ;
		const dreal k_222 = (k_p22 + k_m22) - no2*KS.vx* (k_p22 - k_m22) + KS.vx*KS.vx * k_022 ;

		//Eq G2015(51)
		const dreal C_211 = k_211 - (k_200*k_011 + no2*k_101*k_110)/KS.rho;
		const dreal C_121 = k_121 - (k_020*k_101 + no2*k_110*k_011)/KS.rho;
		const dreal C_112 = k_112 - (k_002*k_110 + no2*k_011*k_101)/KS.rho;

		//Eq G2015(52)
		const dreal C_220 = k_220- (k_020*k_200 + no2*k_110*k_110)/KS.rho;
		const dreal C_022 = k_022- (k_002*k_020 + no2*k_011*k_011)/KS.rho;
		const dreal C_202 = k_202- (k_200*k_002 + no2*k_101*k_101)/KS.rho;

		//Eq G2015(53)
		const dreal C_122 = k_122- (k_020*k_102 + k_002*k_120 + no4*k_011*k_111 + no2*( k_110*k_012 + k_101*k_021))/KS.rho;
		const dreal C_212 = k_212- (k_002*k_210 + k_200*k_012 + no4*k_101*k_111 + no2*( k_011*k_201 + k_110*k_102))/KS.rho;
		const dreal C_221 = k_221- (k_200*k_021 + k_020*k_201 + no4*k_110*k_111 + no2*( k_101*k_120 + k_011*k_210))/KS.rho;
		//Eq G2015(54)
		const dreal C_222 = k_222 - ( no4*k_111*k_111 + k_200*k_022 + k_020*k_202 + k_002*k_220 + no4*(k_011*k_211 + k_101*k_121 + k_110*k_112) + no2*(k_120*k_102 + k_210*k_012 + k_201*k_021))/KS.rho +(no16*k_110*k_101*k_011 + no4*(k_101*k_101*k_020 + k_011*k_011*k_200 + k_110*k_110*k_002)+no2*k_200*k_020*k_002)/KS.rho/KS.rho;
		// gen1nowell.php END

		// relaxation definition
		const dreal omega1 = no1/(no3*KS.lbmViscosity+n1o2); // shear viscosity
		const dreal omega2 = no1;///(no3*KS.lbmViscosity*no2 + n1o2); // bulkViscosity > Viscosity ... test: bulkViscosity = 2/3 shearViscosity
		#ifdef USE_GEIER_CUM_2017
			const dreal lambda3 = (dreal)(0.01); // kolik mam zvolit??? ---> Section 7 @ Geier 2017
			const dreal lambda4 = (dreal)(0.01); // kolik mam zvolit???
			const dreal lambda5 = (dreal)(0.01); // kolik mam zvolit???
			const dreal omega3 = no8*(omega1-no2)*(omega2*(no3*omega1-no1)-no5*omega1)/(no8*(no5-no2*omega1)*omega1 + omega2*(no8+omega1*(no9*omega1-no26)));
			const dreal omega120p102 = omega3 + (no1-omega3)*fabs(C_120+C_102)/(KS.rho*lambda3 + fabs(C_120+C_102)); // limiter
			const dreal omega210p012 = omega3 + (no1-omega3)*fabs(C_120+C_102)/(KS.rho*lambda3 + fabs(C_120+C_102)); // limiter
			const dreal omega201p021 = omega3 + (no1-omega3)*fabs(C_201+C_021)/(KS.rho*lambda3 + fabs(C_201+C_021)); // limiter
			const dreal omega4 = no8*(omega1-no2)*(omega1+omega2*(no3*omega1-no7))/(omega2*(no56-no42*omega1+no8*omega1*omega1)-no8*omega1);
			const dreal omega120m102 = omega4 + (no1-omega4)*fabs(C_120-C_102)/(KS.rho*lambda4 + fabs(C_120-C_102)); // limiter
			const dreal omega210m012 = omega4 + (no1-omega4)*fabs(C_120-C_102)/(KS.rho*lambda4 + fabs(C_120-C_102)); // limiter
			const dreal omega201m021 = omega4 + (no1-omega4)*fabs(C_201-C_021)/(KS.rho*lambda4 + fabs(C_201-C_021)); // limiter
			const dreal omega5 = no24*(omega1-no2)*(no4*omega1*omega1+omega1*omega2*(no18-no13*omega1)+omega2*omega2*(no2+omega1*(no6*omega1-no11)))/(no16*omega1*omega1*(omega1-no6)-no2*omega1*omega2*(no216+no5*omega1*(no9*omega1-no46))+omega2*omega2*(omega1*(no3*omega1-no10)*(no15*omega1-no28)-no48));
			const dreal omega111 = omega5 + (no1-omega5)*fabs(C_111)/(KS.rho*lambda5 + fabs(C_111)); // limiter
			const dreal omega6 = no1;
			const dreal omega7 = no1;
			const dreal omega8 = no1;
			const dreal omega9 = no1;
			const dreal omega10 = no1;
			// extra parameters
			const dreal A=(no4*omega1*omega1 + no2*omega1*omega2*(omega1-no6)+omega2*omega2*(omega1*(no10-no3*omega1)-no4))/(omega1-omega2)/(omega2*(no2+no3*omega1)-no8*omega1);
			const dreal B=(no4*omega1*omega2*(no9*omega1-no16)-no4*omega1*omega1 -no2*omega2*omega2*(no2+no9*omega1*(omega1-no2)))/no3/(omega1-omega2)/(omega2*(no2+no3*omega1)-no8*omega1);
		#else
			const dreal omega3 = no1;
			const dreal omega4 = no1;
			const dreal omega5 = no1;
			const dreal omega6 = no1;
			const dreal omega7 = no1;
			const dreal omega8 = no1;
			const dreal omega9 = no1;
			const dreal omega10 = no1;
			// extra parameters
			const dreal A=0;
			const dreal B=0;
		#endif

		// actual collision step: note: ks = Cs for these indexes
		const dreal Cs_110 = (no1 - omega1)*C_110;
		const dreal Cs_101 = (no1 - omega1)*C_101;
		const dreal Cs_011 = (no1 - omega1)*C_011;
		#ifdef USE_GEIER_CUM_ANTIALIAS
			// derivatives of v: notation taken from Geier's paper 2017 part I: Eq 27-29
//			const dreal Dxu = - omega1/no2/KS.rho * (no2*C_200-C_020-C_002) - omega2/no2/KS.rho*(C_200+C_020+C_002-k_000);
			const dreal Dxu = - omega1/no2/KS.rho * (no2*C_200-C_020-C_002) - omega2/no2/KS.rho*(C_200+C_020+C_002-(no1-KS.rho)); // PE 2019_01_18 poznamka: rho je v tom clanku rho^(2), tj. rho^(2) = 1-rho = 1-k_000
			const dreal Dyv = Dxu + n3o2*omega1/KS.rho *(C_200-C_020);
			const dreal Dzw = Dxu + n3o2*omega1/KS.rho *(C_200-C_002);
			// plus their combination: Eq 30 - 32
			const dreal DxvDyu = -no3*omega1/KS.rho*C_110;
			const dreal DxwDzu = -no3*omega1/KS.rho*C_101;
			const dreal DywDzv = -no3*omega1/KS.rho*C_011;
		#else
			const dreal Dxu = 0;
			const dreal Dyv = 0;
			const dreal Dzw = 0;
			// plus their combination: Eq 30 - 32
			const dreal DxvDyu = 0;
			const dreal DxwDzu = 0;
			const dreal DywDzv = 0;
		#endif
		// Eqs 33-35
		const dreal Eq33RHS = (no1-omega1)*(C_200-C_020) - no3*KS.rho*(no1 - omega1*n1o2)*(KS.vx*KS.vx*Dxu - KS.vy*KS.vy*Dyv);
		const dreal Eq34RHS = (no1-omega1)*(C_200-C_002) - no3*KS.rho*(no1 - omega1*n1o2)*(KS.vx*KS.vx*Dxu - KS.vz*KS.vz*Dzw);
		const dreal Eq35RHS = k_000*omega2 + (no1-omega2)*(C_200+C_020+C_002) - no3*KS.rho*(no1-omega2/no2)*(KS.vx*KS.vx*Dxu + KS.vy*KS.vy*Dyv + KS.vz*KS.vz*Dzw);
		// see rovnice2.mw
		const dreal Cs_200 = n1o3*(Eq33RHS + Eq34RHS + Eq35RHS);
		const dreal Cs_020 = n1o3*(-no2*Eq33RHS + Eq34RHS + Eq35RHS);
		const dreal Cs_002 = n1o3*(Eq33RHS - no2*Eq34RHS + Eq35RHS);

		#ifdef USE_GEIER_CUM_2017
			// Limiter for omegas
			// Eqs 36-41: see vzorce2b.mw
			const dreal Eq117 = (no1 - omega120p102)*(C_120+C_102);
			const dreal Eq118 = (no1 - omega210p012)*(C_210+C_012);
			const dreal Eq119 = (no1 - omega201p021)*(C_201+C_021);
			const dreal Eq120 = (no1 - omega120m102)*(C_120-C_102);
			const dreal Eq121 = (no1 - omega210m012)*(C_210-C_012);
			const dreal Eq122 = (no1 - omega201m021)*(C_201-C_021);

			const dreal Cs_120 = n1o2*(Eq120+Eq117);
			const dreal Cs_102 = n1o2*(-Eq120+Eq117);
			const dreal Cs_210 = n1o2*(Eq121+Eq118);
			const dreal Cs_012 = n1o2*(-Eq121+Eq118);
			const dreal Cs_021 = n1o2*(-Eq122+Eq119);
			const dreal Cs_201 = n1o2*(Eq122+Eq119);
			// Eq 42
			const dreal Cs_111 = (no1-omega111)*C_111;
		#else		
			// Eqs 36-41: see vzorce2.mw
			const dreal Cs_120 = (-C_102-C_120)*omega3*n1o2+(C_102-C_120)*omega4*n1o2+C_120;
			const dreal Cs_102 = (-C_102-C_120)*omega3*n1o2+(-C_102+C_120)*omega4*n1o2+C_102;
			const dreal Cs_210 = (-C_012-C_210)*omega3*n1o2+(C_012-C_210)*omega4*n1o2+C_210;
			const dreal Cs_012 = (-C_012-C_210)*omega3*n1o2+(-C_012+C_210)*omega4*n1o2+C_012;
			const dreal Cs_021 = (-C_021-C_201)*omega3*n1o2+(-C_021+C_201)*omega4*n1o2+C_021;
			const dreal Cs_201 = (-C_021-C_201)*omega3*n1o2+(C_021-C_201)*omega4*n1o2+C_201;
			// Eq 42
			const dreal Cs_111 = (no1-omega5)*C_111;
		#endif		
		// Eqs 43-45
		const dreal Eq43RHS = n2o3*(no1/omega1-n1o2)*omega6*A*KS.rho*(Dxu - no2*Dyv + Dzw) + (no1-omega6)*(C_220-no2*C_202+C_022);
		const dreal Eq44RHS = n2o3*(no1/omega1-n1o2)*omega6*A*KS.rho*(Dxu + Dyv -no2* Dzw) + (no1-omega6)*(C_220+C_202-no2*C_022);
		const dreal Eq45RHS = -n4o3*(no1/omega1-n1o2)*omega7*A*KS.rho*(Dxu + Dyv + Dzw) + (no1-omega7)*(C_220+C_202+C_022);
		// see rovnice2.mw
		const dreal Cs_220 = n1o3*(Eq43RHS+Eq44RHS+Eq45RHS);
		const dreal Cs_202 = n1o3*(-Eq43RHS+Eq45RHS);
		const dreal Cs_022 = n1o3*(-Eq44RHS+Eq45RHS);
		// Eq 46-48
		const dreal Cs_211 = -n1o3*(no1/omega1-n1o2)*omega8*B*KS.rho*DywDzv + (no1-omega8)*C_211;
		const dreal Cs_121 = -n1o3*(no1/omega1-n1o2)*omega8*B*KS.rho*DxwDzu + (no1-omega8)*C_121;
		const dreal Cs_112 = -n1o3*(no1/omega1-n1o2)*omega8*B*KS.rho*DxvDyu + (no1-omega8)*C_112;
		// Eqs 49-52
		const dreal Cs_221 = (no1-omega9)*C_221;
		const dreal Cs_212 = (no1-omega9)*C_212;
		const dreal Cs_122 = (no1-omega9)*C_122;
		const dreal Cs_222 = (no1-omega10)*C_222;

// remark: collision step is the same as in well-CUM

		// 3.4 Backward cumulant transformation
		// gen2nowell.php BEGIN
		//Eq. G2015(81)
		const dreal ks_211 = Cs_211+ (ks_200*ks_011 + no2*ks_101*ks_110)/KS.rho;
		const dreal ks_121 = Cs_121+ (ks_020*ks_101 + no2*ks_110*ks_011)/KS.rho;
		const dreal ks_112 = Cs_112+ (ks_002*ks_110 + no2*ks_011*ks_101)/KS.rho;

		//Eq. G2015(82)
		const dreal ks_220 = Cs_220+ (ks_020*ks_200 + no2*ks_110*ks_110)/KS.rho;
		const dreal ks_022 = Cs_022+ (ks_002*ks_020 + no2*ks_011*ks_011)/KS.rho;
		const dreal ks_202 = Cs_202+ (ks_200*ks_002 + no2*ks_101*ks_101)/KS.rho;

		//Eq. G2015(83)
		const dreal ks_122 = Cs_122+ (ks_020*ks_102 + ks_002*ks_120+ no4*ks_011*ks_111+ no2*(ks_110*ks_012+ks_101*ks_021))/KS.rho;
		const dreal ks_212 = Cs_212+ (ks_002*ks_210 + ks_200*ks_012+ no4*ks_101*ks_111+ no2*(ks_011*ks_201+ks_110*ks_102))/KS.rho;
		const dreal ks_221 = Cs_221+ (ks_200*ks_021 + ks_020*ks_201+ no4*ks_110*ks_111+ no2*(ks_101*ks_120+ks_011*ks_210))/KS.rho;
		// gen2nowell.php END
		
		// Eq. G2015(84)
		const dreal ks_222 = Cs_222 + (no4*ks_111*ks_111 + ks_200*ks_022 + ks_020*ks_202 + ks_002*ks_220 + no4*(ks_011*ks_211 + ks_101*ks_121 + ks_110*ks_112)
					+ no2*(ks_120*ks_102 + ks_210*ks_012 + ks_201*ks_021))/KS.rho
					- (no16*ks_110*ks_101*ks_011 + no4*(ks_101*ks_101*ks_020 + ks_011*ks_011*ks_200 + ks_110*ks_110*ks_002) + no2*ks_200*ks_020*ks_002)/KS.rho/KS.rho;


		// backward central moment transformation
		const dreal ks_000 =  k_000;
		// Geier 2017: forcing scheme ---> zde se musi dat (-1) pro forcing ... ale s nulovou silou a - zde to nefunguje ... 
		const dreal ks_100 =  -k_100;
		const dreal ks_010 =  -k_010;
		const dreal ks_001 =  -k_001;
		
		// gen3nowell.php BEGIN
		//Eq G2015(88)
		const dreal ks_z00 = ks_000*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_100 - ks_200;
		const dreal ks_z01 = ks_001*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_101 - ks_201;
		const dreal ks_z02 = ks_002*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_102 - ks_202;
		const dreal ks_z10 = ks_010*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_110 - ks_210;
		const dreal ks_z11 = ks_011*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_111 - ks_211;
		const dreal ks_z12 = ks_012*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_112 - ks_212;
		const dreal ks_z20 = ks_020*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_120 - ks_220;
		const dreal ks_z21 = ks_021*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_121 - ks_221;
		const dreal ks_z22 = ks_022*(no1-KS.vx*KS.vx) - no2*KS.vx*ks_122 - ks_222;

		//Eq  G2015(89)
		const dreal ks_m00 = (ks_000*(KS.vx*KS.vx - KS.vx) + ks_100*(no2*KS.vx - no1) + ks_200)*n1o2;
		const dreal ks_m01 = (ks_001*(KS.vx*KS.vx - KS.vx) + ks_101*(no2*KS.vx - no1) + ks_201)*n1o2;
		const dreal ks_m02 = (ks_002*(KS.vx*KS.vx - KS.vx) + ks_102*(no2*KS.vx - no1) + ks_202)*n1o2;
		const dreal ks_m10 = (ks_010*(KS.vx*KS.vx - KS.vx) + ks_110*(no2*KS.vx - no1) + ks_210)*n1o2;
		const dreal ks_m11 = (ks_011*(KS.vx*KS.vx - KS.vx) + ks_111*(no2*KS.vx - no1) + ks_211)*n1o2;
		const dreal ks_m12 = (ks_012*(KS.vx*KS.vx - KS.vx) + ks_112*(no2*KS.vx - no1) + ks_212)*n1o2;
		const dreal ks_m20 = (ks_020*(KS.vx*KS.vx - KS.vx) + ks_120*(no2*KS.vx - no1) + ks_220)*n1o2;
		const dreal ks_m21 = (ks_021*(KS.vx*KS.vx - KS.vx) + ks_121*(no2*KS.vx - no1) + ks_221)*n1o2;
		const dreal ks_m22 = (ks_022*(KS.vx*KS.vx - KS.vx) + ks_122*(no2*KS.vx - no1) + ks_222)*n1o2;

		//Eq  G2015(90)
		const dreal ks_p00 = (ks_000*(KS.vx*KS.vx + KS.vx) + ks_100*(no2*KS.vx + no1) + ks_200)*n1o2;
		const dreal ks_p01 = (ks_001*(KS.vx*KS.vx + KS.vx) + ks_101*(no2*KS.vx + no1) + ks_201)*n1o2;
		const dreal ks_p02 = (ks_002*(KS.vx*KS.vx + KS.vx) + ks_102*(no2*KS.vx + no1) + ks_202)*n1o2;
		const dreal ks_p10 = (ks_010*(KS.vx*KS.vx + KS.vx) + ks_110*(no2*KS.vx + no1) + ks_210)*n1o2;
		const dreal ks_p11 = (ks_011*(KS.vx*KS.vx + KS.vx) + ks_111*(no2*KS.vx + no1) + ks_211)*n1o2;
		const dreal ks_p12 = (ks_012*(KS.vx*KS.vx + KS.vx) + ks_112*(no2*KS.vx + no1) + ks_212)*n1o2;
		const dreal ks_p20 = (ks_020*(KS.vx*KS.vx + KS.vx) + ks_120*(no2*KS.vx + no1) + ks_220)*n1o2;
		const dreal ks_p21 = (ks_021*(KS.vx*KS.vx + KS.vx) + ks_121*(no2*KS.vx + no1) + ks_221)*n1o2;
		const dreal ks_p22 = (ks_022*(KS.vx*KS.vx + KS.vx) + ks_122*(no2*KS.vx + no1) + ks_222)*n1o2;

		//Eq G2015(91)
		const dreal ks_mz0 = ks_m00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m10 - ks_m20;
		const dreal ks_mz1 = ks_m01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m11 - ks_m21;
		const dreal ks_mz2 = ks_m02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_m12 - ks_m22;
		const dreal ks_zz0 = ks_z00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z10 - ks_z20;
		const dreal ks_zz1 = ks_z01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z11 - ks_z21;
		const dreal ks_zz2 = ks_z02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_z12 - ks_z22;
		const dreal ks_pz0 = ks_p00*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p10 - ks_p20;
		const dreal ks_pz1 = ks_p01*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p11 - ks_p21;
		const dreal ks_pz2 = ks_p02*(no1 - KS.vy*KS.vy) - no2*KS.vy*ks_p12 - ks_p22;

		//Eq  G2015(92)
		const dreal ks_mm0 = (ks_m00*(KS.vy*KS.vy - KS.vy) + ks_m10*(no2*KS.vy - no1) + ks_m20)*n1o2;
		const dreal ks_mm1 = (ks_m01*(KS.vy*KS.vy - KS.vy) + ks_m11*(no2*KS.vy - no1) + ks_m21)*n1o2;
		const dreal ks_mm2 = (ks_m02*(KS.vy*KS.vy - KS.vy) + ks_m12*(no2*KS.vy - no1) + ks_m22)*n1o2;
		const dreal ks_zm0 = (ks_z00*(KS.vy*KS.vy - KS.vy) + ks_z10*(no2*KS.vy - no1) + ks_z20)*n1o2;
		const dreal ks_zm1 = (ks_z01*(KS.vy*KS.vy - KS.vy) + ks_z11*(no2*KS.vy - no1) + ks_z21)*n1o2;
		const dreal ks_zm2 = (ks_z02*(KS.vy*KS.vy - KS.vy) + ks_z12*(no2*KS.vy - no1) + ks_z22)*n1o2;
		const dreal ks_pm0 = (ks_p00*(KS.vy*KS.vy - KS.vy) + ks_p10*(no2*KS.vy - no1) + ks_p20)*n1o2;
		const dreal ks_pm1 = (ks_p01*(KS.vy*KS.vy - KS.vy) + ks_p11*(no2*KS.vy - no1) + ks_p21)*n1o2;
		const dreal ks_pm2 = (ks_p02*(KS.vy*KS.vy - KS.vy) + ks_p12*(no2*KS.vy - no1) + ks_p22)*n1o2;

		//Eq G2015(93)
		const dreal ks_mp0 = (ks_m00*(KS.vy*KS.vy + KS.vy) + ks_m10*(no2*KS.vy + no1) + ks_m20)*n1o2;
		const dreal ks_mp1 = (ks_m01*(KS.vy*KS.vy + KS.vy) + ks_m11*(no2*KS.vy + no1) + ks_m21)*n1o2;
		const dreal ks_mp2 = (ks_m02*(KS.vy*KS.vy + KS.vy) + ks_m12*(no2*KS.vy + no1) + ks_m22)*n1o2;
		const dreal ks_zp0 = (ks_z00*(KS.vy*KS.vy + KS.vy) + ks_z10*(no2*KS.vy + no1) + ks_z20)*n1o2;
		const dreal ks_zp1 = (ks_z01*(KS.vy*KS.vy + KS.vy) + ks_z11*(no2*KS.vy + no1) + ks_z21)*n1o2;
		const dreal ks_zp2 = (ks_z02*(KS.vy*KS.vy + KS.vy) + ks_z12*(no2*KS.vy + no1) + ks_z22)*n1o2;
		const dreal ks_pp0 = (ks_p00*(KS.vy*KS.vy + KS.vy) + ks_p10*(no2*KS.vy + no1) + ks_p20)*n1o2;
		const dreal ks_pp1 = (ks_p01*(KS.vy*KS.vy + KS.vy) + ks_p11*(no2*KS.vy + no1) + ks_p21)*n1o2;
		const dreal ks_pp2 = (ks_p02*(KS.vy*KS.vy + KS.vy) + ks_p12*(no2*KS.vy + no1) + ks_p22)*n1o2;

		//Eq G2015(94)
		KS.f[mmz] = ks_mm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mm1 - ks_mm2;
		KS.f[mzz] = ks_mz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mz1 - ks_mz2;
		KS.f[mpz] = ks_mp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_mp1 - ks_mp2;
		KS.f[zmz] = ks_zm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zm1 - ks_zm2;
		KS.f[zzz] = ks_zz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zz1 - ks_zz2;
		KS.f[zpz] = ks_zp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_zp1 - ks_zp2;
		KS.f[pmz] = ks_pm0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pm1 - ks_pm2;
		KS.f[pzz] = ks_pz0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pz1 - ks_pz2;
		KS.f[ppz] = ks_pp0*(no1-KS.vz*KS.vz) - no2*KS.vz*ks_pp1 - ks_pp2;

		//Eq  G2015(95)
		KS.f[mmm] = (ks_mm0*(KS.vz*KS.vz-KS.vz) + ks_mm1*(no2*KS.vz - no1) + ks_mm2)*n1o2;
		KS.f[mzm] = (ks_mz0*(KS.vz*KS.vz-KS.vz) + ks_mz1*(no2*KS.vz - no1) + ks_mz2)*n1o2;
		KS.f[mpm] = (ks_mp0*(KS.vz*KS.vz-KS.vz) + ks_mp1*(no2*KS.vz - no1) + ks_mp2)*n1o2;
		KS.f[zmm] = (ks_zm0*(KS.vz*KS.vz-KS.vz) + ks_zm1*(no2*KS.vz - no1) + ks_zm2)*n1o2;
		KS.f[zzm] = (ks_zz0*(KS.vz*KS.vz-KS.vz) + ks_zz1*(no2*KS.vz - no1) + ks_zz2)*n1o2;
		KS.f[zpm] = (ks_zp0*(KS.vz*KS.vz-KS.vz) + ks_zp1*(no2*KS.vz - no1) + ks_zp2)*n1o2;
		KS.f[pmm] = (ks_pm0*(KS.vz*KS.vz-KS.vz) + ks_pm1*(no2*KS.vz - no1) + ks_pm2)*n1o2;
		KS.f[pzm] = (ks_pz0*(KS.vz*KS.vz-KS.vz) + ks_pz1*(no2*KS.vz - no1) + ks_pz2)*n1o2;
		KS.f[ppm] = (ks_pp0*(KS.vz*KS.vz-KS.vz) + ks_pp1*(no2*KS.vz - no1) + ks_pp2)*n1o2;

		//Eq  G2015(96)
		KS.f[mmp] = (ks_mm0*(KS.vz*KS.vz+KS.vz) + ks_mm1*(no2*KS.vz + no1) + ks_mm2)*n1o2;
		KS.f[mzp] = (ks_mz0*(KS.vz*KS.vz+KS.vz) + ks_mz1*(no2*KS.vz + no1) + ks_mz2)*n1o2;
		KS.f[mpp] = (ks_mp0*(KS.vz*KS.vz+KS.vz) + ks_mp1*(no2*KS.vz + no1) + ks_mp2)*n1o2;
		KS.f[zmp] = (ks_zm0*(KS.vz*KS.vz+KS.vz) + ks_zm1*(no2*KS.vz + no1) + ks_zm2)*n1o2;
		KS.f[zzp] = (ks_zz0*(KS.vz*KS.vz+KS.vz) + ks_zz1*(no2*KS.vz + no1) + ks_zz2)*n1o2;
		KS.f[zpp] = (ks_zp0*(KS.vz*KS.vz+KS.vz) + ks_zp1*(no2*KS.vz + no1) + ks_zp2)*n1o2;
		KS.f[pmp] = (ks_pm0*(KS.vz*KS.vz+KS.vz) + ks_pm1*(no2*KS.vz + no1) + ks_pm2)*n1o2;
		KS.f[pzp] = (ks_pz0*(KS.vz*KS.vz+KS.vz) + ks_pz1*(no2*KS.vz + no1) + ks_pz2)*n1o2;
		KS.f[ppp] = (ks_pp0*(KS.vz*KS.vz+KS.vz) + ks_pp1*(no2*KS.vz + no1) + ks_pp2)*n1o2;
		// gen3nowell.php END
		#undef C_000
		#undef C_011
		#undef C_010
		#undef C_100
		#undef C_101
		#undef C_110
		#undef C_111
		#undef C_002
		#undef C_020
		#undef C_200
		#undef C_012
		#undef C_021
		#undef C_102
		#undef C_201
		#undef C_120
		#undef C_210
		// backward translation
		#undef ks_000
		#undef ks_001
		#undef ks_011
		#undef ks_010
		#undef ks_100
		#undef ks_101
		#undef ks_110
		#undef ks_111
		#undef ks_002
		#undef ks_020
		#undef ks_200
		#undef ks_012
		#undef ks_021
		#undef ks_102
		#undef ks_201
		#undef ks_120
		#undef ks_210
		// backward translation
		#undef Cs_001
		#undef Cs_010
		#undef Cs_100
	}
};
