#include "lbm_common.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_FCLBM : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "FCLBM"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
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
		const dreal k_mm1 = (KS.f[mmp] - KS.f[mmm]) - KS.vz * ( k_mm0 );
		const dreal k_mz1 = (KS.f[mzp] - KS.f[mzm]) - KS.vz * ( k_mz0 );
		const dreal k_mp1 = (KS.f[mpp] - KS.f[mpm]) - KS.vz * ( k_mp0 );
		const dreal k_zm1 = (KS.f[zmp] - KS.f[zmm]) - KS.vz * ( k_zm0 );
		const dreal k_zz1 = (KS.f[zzp] - KS.f[zzm]) - KS.vz * ( k_zz0 );
		const dreal k_zp1 = (KS.f[zpp] - KS.f[zpm]) - KS.vz * ( k_zp0 );
		const dreal k_pm1 = (KS.f[pmp] - KS.f[pmm]) - KS.vz * ( k_pm0 );
		const dreal k_pz1 = (KS.f[pzp] - KS.f[pzm]) - KS.vz * ( k_pz0 );
		const dreal k_pp1 = (KS.f[ppp] - KS.f[ppm]) - KS.vz * ( k_pp0 );

		//Eq 8
		const dreal k_mm2 = (KS.f[mmp] + KS.f[mmm]) - no2*KS.vz*(KS.f[mmp] - KS.f[mmm]) + KS.vz*KS.vz*( k_mm0 );
		const dreal k_mz2 = (KS.f[mzp] + KS.f[mzm]) - no2*KS.vz*(KS.f[mzp] - KS.f[mzm]) + KS.vz*KS.vz*( k_mz0 );
		const dreal k_mp2 = (KS.f[mpp] + KS.f[mpm]) - no2*KS.vz*(KS.f[mpp] - KS.f[mpm]) + KS.vz*KS.vz*( k_mp0 );
		const dreal k_zm2 = (KS.f[zmp] + KS.f[zmm]) - no2*KS.vz*(KS.f[zmp] - KS.f[zmm]) + KS.vz*KS.vz*( k_zm0 );
		const dreal k_zz2 = (KS.f[zzp] + KS.f[zzm]) - no2*KS.vz*(KS.f[zzp] - KS.f[zzm]) + KS.vz*KS.vz*( k_zz0 );
		const dreal k_zp2 = (KS.f[zpp] + KS.f[zpm]) - no2*KS.vz*(KS.f[zpp] - KS.f[zpm]) + KS.vz*KS.vz*( k_zp0 );
		const dreal k_pm2 = (KS.f[pmp] + KS.f[pmm]) - no2*KS.vz*(KS.f[pmp] - KS.f[pmm]) + KS.vz*KS.vz*( k_pm0 );
		const dreal k_pz2 = (KS.f[pzp] + KS.f[pzm]) - no2*KS.vz*(KS.f[pzp] - KS.f[pzm]) + KS.vz*KS.vz*( k_pz0 );
		const dreal k_pp2 = (KS.f[ppp] + KS.f[ppm]) - no2*KS.vz*(KS.f[ppp] - KS.f[ppm]) + KS.vz*KS.vz*( k_pp0 );

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
		const dreal k_m10 = (k_mp0 - k_mm0) - KS.vy * (  k_m00 );
		const dreal k_z10 = (k_zp0 - k_zm0) - KS.vy * (  k_z00 );
		const dreal k_p10 = (k_pp0 - k_pm0) - KS.vy * (  k_p00 );
		const dreal k_m11 = (k_mp1 - k_mm1) - KS.vy * (  k_m01 );
		const dreal k_z11 = (k_zp1 - k_zm1) - KS.vy * (  k_z01 );
		const dreal k_p11 = (k_pp1 - k_pm1) - KS.vy * (  k_p01 );
		const dreal k_m12 = (k_mp2 - k_mm2) - KS.vy * (  k_m02 );
		const dreal k_z12 = (k_zp2 - k_zm2) - KS.vy * (  k_z02 );
		const dreal k_p12 = (k_pp2 - k_pm2) - KS.vy * (  k_p02 );

		//Eq 11
		const dreal k_m20 = (k_mp0 + k_mm0) - no2*KS.vy* (k_mp0 - k_mm0) + KS.vy*KS.vy * (  k_m00 );
		const dreal k_z20 = (k_zp0 + k_zm0) - no2*KS.vy* (k_zp0 - k_zm0) + KS.vy*KS.vy * (  k_z00 );
		const dreal k_p20 = (k_pp0 + k_pm0) - no2*KS.vy* (k_pp0 - k_pm0) + KS.vy*KS.vy * (  k_p00  );
		const dreal k_m21 = (k_mp1 + k_mm1) - no2*KS.vy* (k_mp1 - k_mm1) + KS.vy*KS.vy * (  k_m01 );
		const dreal k_z21 = (k_zp1 + k_zm1) - no2*KS.vy* (k_zp1 - k_zm1) + KS.vy*KS.vy * (  k_z01 );
		const dreal k_p21 = (k_pp1 + k_pm1) - no2*KS.vy* (k_pp1 - k_pm1) + KS.vy*KS.vy * (  k_p01 );
		const dreal k_m22 = (k_mp2 + k_mm2) - no2*KS.vy* (k_mp2 - k_mm2) + KS.vy*KS.vy * (  k_m02  );
		const dreal k_z22 = (k_zp2 + k_zm2) - no2*KS.vy* (k_zp2 - k_zm2) + KS.vy*KS.vy * (  k_z02 );
		const dreal k_p22 = (k_pp2 + k_pm2) - no2*KS.vy* (k_pp2 - k_pm2) + KS.vy*KS.vy * (  k_p02 );

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
		
		// bacha: k_000 neni rho, ale je to delta_rho

		//Eq 13
		const dreal k_100 = (k_p00 - k_m00) - KS.vx * (  k_000 );
		const dreal k_101 = (k_p01 - k_m01) - KS.vx * (  k_001 );
		const dreal k_102 = (k_p02 - k_m02) - KS.vx * (  k_002 );
		const dreal k_110 = (k_p10 - k_m10) - KS.vx * (  k_010 );
		const dreal k_111 = (k_p11 - k_m11) - KS.vx * (  k_011 );
		const dreal k_112 = (k_p12 - k_m12) - KS.vx * (  k_012 );
		const dreal k_120 = (k_p20 - k_m20) - KS.vx * (  k_020 );
		const dreal k_121 = (k_p21 - k_m21) - KS.vx * (  k_021 );
		const dreal k_122 = (k_p22 - k_m22) - KS.vx * (  k_022 );

		//Eq 14
		const dreal k_200 = (k_p00 + k_m00) - no2*KS.vx* (k_p00 - k_m00) + KS.vx*KS.vx * (  k_000 );
		const dreal k_201 = (k_p01 + k_m01) - no2*KS.vx* (k_p01 - k_m01) + KS.vx*KS.vx * (  k_001 );
		const dreal k_202 = (k_p02 + k_m02) - no2*KS.vx* (k_p02 - k_m02) + KS.vx*KS.vx * (  k_002 );
		const dreal k_210 = (k_p10 + k_m10) - no2*KS.vx* (k_p10 - k_m10) + KS.vx*KS.vx * (  k_010 );
		const dreal k_211 = (k_p11 + k_m11) - no2*KS.vx* (k_p11 - k_m11) + KS.vx*KS.vx * (  k_011 );
		const dreal k_212 = (k_p12 + k_m12) - no2*KS.vx* (k_p12 - k_m12) + KS.vx*KS.vx * (  k_012 );
		const dreal k_220 = (k_p20 + k_m20) - no2*KS.vx* (k_p20 - k_m20) + KS.vx*KS.vx * (  k_020 );
		const dreal k_221 = (k_p21 + k_m21) - no2*KS.vx* (k_p21 - k_m21) + KS.vx*KS.vx * (  k_021 );
		const dreal k_222 = (k_p22 + k_m22) - no2*KS.vx* (k_p22 - k_m22) + KS.vx*KS.vx * (  k_022 );
		// tak, tedka mam centralni momenty
		// gen1.php END

		// relaxation definition
		dreal omega1 = no1/(no3*KS.lbmViscosity+n1o2);
		dreal omega2 = no1;///(no3*KS.lbmViscosity+n1o2); // enhanced bulk viscosity
		dreal omega3 = no1;
		dreal omega4 = no1;
		dreal omega5 = no1;
		dreal omega6 = no1;
		dreal omega7 = no1;
		dreal omega8 = no1;
		dreal omega9 = no1;
		dreal omega10 = no1;

		dreal omega [23];
		omega[0] = omega1;
		omega[1] = omega1;
		omega[2] = omega1;
		omega[3] = omega1;
		omega[4] = omega1;

		for(int i = 5; i <= 22; i++)
		    omega[i] = (dreal)1.;

		
		dreal g_4, g_5, g_6, g_7, g_8, g_9, g_10, g_11, g_12, g_13, g_14, g_15, g_16, g_17, g_18, g_19, g_20, g_21, g_22, g_23, g_24, g_25, g_26;
		
		g_4 = n1o12*omega[0]*(-k_110);
		g_5 = n1o12*omega[1]*(-k_101);
		g_6 = n1o12*omega[2]*(-k_011);
		
		g_7 = n1o12*omega[3]*(-k_200 + k_020);
		g_8 = n1o36*omega[4]*(-1.*(k_200 + k_020 - (dreal)2.0*k_002));
		g_9 = n1o18*omega[5]*(-1.*(k_200 + k_020 + k_002) + KS.rho);
		
		g_10 = n1o24*omega[6]*( -k_102 - k_120) + KS.vy*g_4 + g_5*KS.vz - n1o4*KS.vx*g_7 - n1o4*KS.vx*g_8 + n1o2*KS.vx*g_9;
		g_11 = n1o24*omega[7]*(-k_012 - k_210) + g_4*KS.vx + KS.vz*g_6 + n1o4*g_7*KS.vy - n1o4*KS.vy*g_8 + n1o2*KS.vy*g_9;
		g_12 = n1o24*omega[8]*(-k_201 - k_021) + KS.vx*g_5 + KS.vy*g_6 + n1o2*KS.vz*(g_8 + g_9);
		
		g_13 = n1o8*omega[9]*( -k_120 + k_102) + no3*(KS.vy*g_4 - KS.vz*g_5) + n3o4*KS.vx*(-g_7 + no3*g_8);
		g_14 = n1o8*omega[10]*(-k_210 + k_012) + no3*(KS.vx*g_4 - KS.vz*g_6) + n3o4*KS.vy*(g_7 + no3*g_8);
		g_15 = n1o8*omega[11]*(-k_201 + k_021) + no3*(KS.vx*g_5 - KS.vy*g_6) + n3o2*KS.vz*g_7;
		g_16 = n1o8*omega[12]*(-k_111) + n3o2*(KS.vz*g_4 + KS.vy*g_5 + KS.vx*g_6);
		
		
		//post collision moments:
		
		dreal kp_200, kp_020, kp_002;
		
		kp_200 = k_200 + (dreal)6.*(g_7 + g_8 + g_9);
		kp_020 = k_020 + (dreal)6.*(-g_7 + g_8 + g_9);
		kp_002 = k_002 - (dreal)12.*g_8 + (dreal)6.*g_9;
		
		
		g_17 = n1o12*omega[13]*(-k_220 - k_202 - k_022 + (kp_200*kp_020 + kp_200*kp_002 + kp_020*kp_002)) - no4*(KS.vx*KS.vy*g_4 + KS.vx*KS.vz*g_5 + KS.vy*KS.vz*g_6) 
			+ n1o2*(KS.vx*KS.vx - KS.vy*KS.vy)*g_7 + n1o2*(KS.vx*KS.vx + KS.vy*KS.vy - no2*KS.vz*KS.vz)*g_8 + (-KS.vx*KS.vx - KS.vy*KS.vy - KS.vz*KS.vz - (dreal)2.)*g_9
			+ no4*(KS.vx*g_10 + KS.vy*g_11 + KS.vz*g_12);
			
		g_18 = n1o24*omega[14]*(-k_220 - k_202 + no2*k_022 + (kp_200*kp_020 + kp_200*kp_002 - no2*kp_020*kp_002)) - no2*KS.vx*KS.vy*g_4
			- no2*KS.vx*KS.vz*g_5 + no4*KS.vy*KS.vz*g_6 + n1o4*(KS.vx*KS.vx - KS.vy*KS.vy - no3*KS.vz*KS.vz - no2)*g_7 
			+ n1o4*(KS.vx*KS.vx - no5*KS.vy*KS.vy + KS.vz*KS.vz - no2)*g_8  + n1o4*(-no2*KS.vx*KS.vx + KS.vy*KS.vy + KS.vz*KS.vz)*g_9
			+ no2*KS.vx*g_10 - KS.vy*g_11 - KS.vz*g_12 + KS.vy*g_14 + KS.vz*g_15;
			
		g_19 = n1o8*omega[15]*(-k_220 + k_202 + (kp_200*kp_020 - kp_200*kp_002)) - no6*(KS.vx*KS.vy*g_4 - KS.vx*KS.vz*g_5) 
			+ n1o4*(no3*KS.vx*KS.vx - no3*KS.vy*KS.vy + no3*KS.vz*KS.vz + no2)*g_7
			+ n1o4*(-no9*KS.vx*KS.vx - no3*KS.vy*KS.vy + no3*KS.vz*KS.vz - no6)*g_8 + n3o4*(-KS.vy*KS.vy + KS.vz*KS.vz)*g_9
			+ no3*(KS.vy*g_11 - KS.vz*g_12) + no2*KS.vx*g_13 + KS.vy*g_14 - KS.vz*g_15;
			
		
		g_20 = n1o8*omega[16]*(-k_211) - no3*(KS.vx*KS.vz*g_4 + KS.vx*KS.vy*g_5) - (n3o2*KS.vx*KS.vx + no1)*g_6 
			- n3o4*(KS.vy*KS.vz*g_7 + KS.vy*KS.vz*g_8 + KS.vy*KS.vz*g_9) + n3o2*(KS.vz*g_11 + KS.vy*g_12) + n1o2*KS.vz*g_14
			+ n1o2*KS.vy*g_15 + (dreal)2.*KS.vx*g_16;
			
			
		g_21 = n1o8*omega[17]*(-k_121) - no3*KS.vy*KS.vz*g_4 - (n3o2*KS.vy*KS.vy + no1)*g_5 - no3*KS.vx*KS.vy*g_6
			-n3o4*(-KS.vx*KS.vz*g_7 + KS.vx*KS.vz*g_8 + KS.vx*KS.vz*g_9) + n3o2*(KS.vz*g_10 + KS.vx*g_12) + n1o2*KS.vz*g_13 - n1o2*KS.vx*g_15
			+ no2*KS.vy*g_16;
			
			
		g_22 = n1o8*omega[18]*(-k_112) - (n3o2*KS.vz*KS.vz + no1)*g_4 - no3*(KS.vy*KS.vz*g_5 + KS.vx*KS.vz*g_6) + n3o2*KS.vx*KS.vy*g_8
			- n3o4*(KS.vx*KS.vy*g_9 - no2*KS.vy*g_10 - no2*KS.vx*g_11) - n1o2*(KS.vy*g_13 + KS.vx*g_14) + no2*KS.vz*g_16;
			
			
		g_23 = 	n1o8*omega[19]*(-k_122) + (no3*KS.vy*KS.vz*KS.vz + no2*KS.vy)*g_4 + (no3*KS.vy*KS.vy*KS.vz + no2*KS.vz)*g_5 + no6*KS.vx*KS.vy*KS.vz*g_6
			- (n3o4*KS.vx*KS.vz*KS.vz + n1o2*KS.vx)*g_7 - n1o4*(no6*KS.vx*KS.vy*KS.vy - no3*KS.vx*KS.vz*KS.vz + no2*KS.vx)*g_8
			+ (n3o4*KS.vx*KS.vy*KS.vy + n3o4*KS.vx*KS.vz*KS.vz + KS.vx)*g_9 - (n3o2*KS.vy*KS.vy + n3o2*KS.vz*KS.vz + no2)*g_10
			 - no3*(KS.vx*KS.vy*g_11 + KS.vx*KS.vz*g_12) + (n1o2*KS.vy*KS.vy - n1o2*KS.vz*KS.vz)*g_13 + KS.vx*KS.vy*g_14 + KS.vx*KS.vz*g_15 - no4*KS.vy*KS.vz*g_16
			 + n1o2*KS.vx*g_17 - KS.vx*g_18 + no2*(KS.vz*g_21 + KS.vy*g_22);
			 
			 
		g_24 = 	 n1o8*omega[20]*(-k_212) + (no3*KS.vx*KS.vz*KS.vz + no2*KS.vx)*g_4 + no6*KS.vx*KS.vy*KS.vz*g_5 + (no3*KS.vx*KS.vx*KS.vz + no2*KS.vz)*g_6
			+ (n3o4*KS.vy*KS.vz*KS.vz + n1o2*KS.vy)*g_7 + (-n3o2*KS.vx*KS.vx*KS.vy + n3o4*KS.vy*KS.vz*KS.vz - n1o2*KS.vy)*g_8
			+ (n3o4*(KS.vx*KS.vx*KS.vy + KS.vy*KS.vz*KS.vz) + KS.vy)*g_9 - no3*KS.vx*KS.vy*g_10 - (n3o2*(KS.vx*KS.vx + KS.vz*KS.vz) + no2)*g_11 
			- no3*KS.vy*KS.vz*g_12 + KS.vx*KS.vy*g_13 + n1o2*(KS.vx*KS.vx - KS.vz*KS.vz)*g_14 - KS.vy*KS.vz*g_15 - no4*KS.vx*KS.vz*g_16 + n1o2*KS.vy*g_17 + n1o2*KS.vy*g_18
			- n1o2*KS.vy*g_19 + no2*KS.vz*g_20 + no2*KS.vx*g_22;
			
			
		g_25 = n1o8*omega[21]*(-k_221) + no6*KS.vx*KS.vy*KS.vz*g_4 + (no3*KS.vx*KS.vy*KS.vy + no2*KS.vx)*g_5 + (no3*KS.vx*KS.vx*KS.vy + no2*KS.vy)*g_6
			+ n3o4*(-KS.vx*KS.vx*KS.vz + KS.vy*KS.vy*KS.vz)*g_7 + (n3o4*KS.vx*KS.vx*KS.vz + n3o4*KS.vy*KS.vy*KS.vz + KS.vz)*g_8
			+ (n3o4*KS.vx*KS.vx*KS.vz + n3o4*KS.vy*KS.vy*KS.vz + KS.vz)*g_9 - no3*(KS.vx*KS.vz*g_10 + KS.vy*KS.vz*g_11)
			-(n3o2*KS.vx*KS.vx + n3o2*KS.vy*KS.vy + no2)*g_12 - KS.vx*KS.vz*g_13 - KS.vy*KS.vz*g_14 + (n1o2*KS.vx*KS.vx - n1o2*KS.vy*KS.vy)*g_15
			- no4*KS.vx*KS.vy*g_16 + n1o2*(KS.vz*g_17 + KS.vz*g_18 + KS.vz*g_19) + no2*(KS.vy*g_20 + KS.vx*g_21);
			
			
		g_26 = n1o8*omega[22]*(kp_200*kp_020*kp_002 - k_222) - (no4*KS.vx*KS.vy + no6*KS.vx*KS.vy*KS.vz*KS.vz)*g_4 - (no6*KS.vx*KS.vy*KS.vy*KS.vz + no4*KS.vx*KS.vz)*g_5
			- (no6*KS.vx*KS.vx*KS.vy*KS.vz + no4*KS.vy*KS.vz)*g_6 + (n1o2*(KS.vx*KS.vx - KS.vy*KS.vy) + n3o4*(KS.vx*KS.vx*KS.vz*KS.vz - KS.vy*KS.vy*KS.vz*KS.vz))*g_7
			+ (n1o2*KS.vx*KS.vx + n1o2*KS.vy*KS.vy - KS.vz*KS.vz + n3o4*(no2*KS.vx*KS.vx*KS.vy*KS.vy - KS.vx*KS.vx*KS.vz*KS.vz - KS.vy*KS.vy*KS.vz*KS.vz))*g_8
			+ (-KS.vx*KS.vx - KS.vy*KS.vy -KS.vz*KS.vz - n3o4*(KS.vx*KS.vx*KS.vy*KS.vy + KS.vx*KS.vx*KS.vz*KS.vz + KS.vy*KS.vy*KS.vz*KS.vz) - no1)*g_9
			+ (no3*(KS.vx*KS.vy*KS.vy + KS.vx*KS.vz*KS.vz) + no4*KS.vx)*g_10 + (no3*(KS.vx*KS.vx*KS.vy + KS.vy*KS.vz*KS.vz) + no4*KS.vy)*g_11
			+ (no3*(KS.vx*KS.vx*KS.vz + KS.vy*KS.vy*KS.vz) + no4*KS.vz)*g_12 + (KS.vx*KS.vz*KS.vz - KS.vx*KS.vy*KS.vy)*g_13 + (KS.vy*KS.vz*KS.vz - KS.vx*KS.vx*KS.vy)*g_14
			+ (KS.vy*KS.vy*KS.vz - KS.vx*KS.vx*KS.vz)*g_15 + no8*KS.vx*KS.vy*KS.vz*g_16 + (n1o2*(-KS.vx*KS.vx - KS.vy*KS.vy - KS.vz*KS.vz) - no1)*g_17
			+ (KS.vx*KS.vx - n1o2*(KS.vy*KS.vy + KS.vz*KS.vz))*g_18 + n1o2*(KS.vy*KS.vy - KS.vz*KS.vz)*g_19 - no4*(KS.vy*KS.vz*g_20 + KS.vx*KS.vz*g_21 + KS.vx*KS.vy*g_22)
			+ no2*(KS.vx*g_23 + KS.vy*g_24 + KS.vz*g_25);
		
//		if(isnan(g_26))
//		printf("g_26 = %1.10f   ", g_26);
		
		
		
		dreal m_0, m_1, m_2, m_3, m_4, m_5, m_6, m_7, m_8, m_9, m_10, m_11, m_12, m_13, m_14, m_15, m_16, m_17, m_18, m_19, m_20, m_21, m_22, m_23, m_24, m_25, m_26;
        
        m_0 = 0;
        m_1 = KS.fx;
        m_2 = KS.fy;
        m_3 = KS.fz;
        m_4 = (KS.fx*KS.vy + KS.fy*KS.vx);
        m_5 = (KS.fx*KS.vz + KS.fz*KS.vx);
        m_6 = (KS.fy*KS.vz + KS.fz*KS.vy);
        m_7 = no2*(KS.fx*KS.vx - KS.fy*KS.vy);
        m_8 = no2*(KS.fx*KS.vx + KS.fy*KS.vy - no2*KS.fz*KS.vz);
        m_9 = no2*(KS.fx*KS.vx + KS.fy*KS.vy + KS.fz*KS.vz);
        m_10 = (no3*KS.vy*KS.vy + no3*KS.vz*KS.vz - no4)*KS.fx + no6*KS.vx*KS.vy*KS.fy + no6*KS.vx*KS.vz*KS.fz;
        m_11 = no6*KS.vx*KS.vy*KS.fx + (no3*KS.vx*KS.vx + no3*KS.vz*KS.vz - no4)*KS.fy + no6*KS.vz*KS.vy*KS.fz;
        m_12 = no6*KS.vx*KS.vz*KS.fx + no6*KS.vz*KS.vy*KS.fy + (no3*KS.vx*KS.vx + no3*KS.vy*KS.vy - no4)*KS.fz;
        m_13 = (KS.vy*KS.vy - KS.vz*KS.vz)*KS.fx + no2*KS.vx*KS.vy*KS.fy - no2*KS.vx*KS.vz*KS.fz;
        m_14 = no2*KS.vx*KS.vy*KS.fx + (KS.vx*KS.vx - KS.vz*KS.vz)*KS.fy - no2*KS.vz*KS.vy*KS.fz;
        m_15 = no2*KS.vx*KS.vz*KS.fx - no2*KS.vz*KS.vy*KS.fy + (KS.vx*KS.vx - KS.vy*KS.vy)*KS.fz;
        m_16 = KS.fx*KS.vy*KS.vz + KS.fy*KS.vx*KS.vz + KS.fz*KS.vx*KS.vy;
        m_17 = (no6*KS.vy*KS.vy + no6*KS.vz*KS.vz - no8)*KS.vx*KS.fx + (no6*KS.vx*KS.vx*KS.vy + no6*KS.vy*KS.vz*KS.vz - no8*KS.vy)*KS.fy + (no6*KS.vx*KS.vx*KS.vz + 6*KS.vy*KS.vy*KS.vz - no8*KS.vz)*KS.fz;
        m_18 = (no6*KS.vy*KS.vy + no6*KS.vz*KS.vz - no8)*KS.vx*KS.fx + (no6*KS.vx*KS.vx*KS.vy - no12*KS.vy*KS.vz*KS.vz + no4*KS.vy)*KS.fy + (no6*KS.vx*KS.vx*KS.vz - no12*KS.vy*KS.vy*KS.vz + no4*KS.vz)*KS.fz;
        m_19 = (no6*KS.vy*KS.vy - no6*KS.vz*KS.vz)*KS.vx*KS.fx + (no6*KS.vx*KS.vx*KS.vy - no4*KS.vy)*KS.fy + (-no6*KS.vx*KS.vx*KS.vz + no4*KS.vz)*KS.fz;
        m_20 = no6*KS.vx*KS.vy*KS.vz*KS.fx + (no3*KS.vx*KS.vx*KS.vz - no2*KS.vz)*KS.fy + (no3*KS.vx*KS.vx*KS.vy - no2*KS.vy)*KS.fz;
        m_21 = (no3*KS.vy*KS.vy*KS.vz - no2*KS.vz)*KS.fx + no6*KS.vx*KS.vy*KS.vz*KS.fy + (no3*KS.vx*KS.vy*KS.vy - no2*KS.vx)*KS.fz;
        m_22 = (no3*KS.vy*KS.vz*KS.vz - no2*KS.vy)*KS.fx + (no3*KS.vx*KS.vz*KS.vz - no2*KS.vx)*KS.fy + no6*KS.vx*KS.vy*KS.vz*KS.fz;
        m_23 = ((no9*KS.vz*KS.vz - no6)*KS.vy*KS.vy - no6*KS.vz*KS.vz + no4)*KS.fx + (no18*KS.vz*KS.vz - no12)*KS.vy*KS.vx*KS.fy + no6*KS.vx*KS.vz*(no3*KS.vy*KS.vy - no2)*KS.fz;
        m_24 = (no18*KS.vz*KS.vz - no12)*KS.vy*KS.vx*KS.fx + ((no9*KS.vz*KS.vz - no6)*KS.vx*KS.vx - no6*KS.vz*KS.vz + no4)*KS.fy + no6*KS.vz*KS.vy*(no3*KS.vx*KS.vx - no2)*KS.fz;
        m_25 = no6*KS.vx*KS.vz*(no3*KS.vy*KS.vy - no2)*KS.fx + no6*KS.vz*KS.vy*(no3*KS.vx*KS.vx - no2)*KS.fy + ((no9*KS.vy*KS.vy - no6)*KS.vx*KS.vx - no6*KS.vy*KS.vy + no4)*KS.fz;
        m_26 = (no6*(no3*KS.vz*KS.vz - no2))*(no3*KS.vy*KS.vy - no2)*KS.vx*KS.fx + (no6*(no3*KS.vz*KS.vz - no2))*(no3*KS.vx*KS.vx - no2)*KS.vy*KS.fy + no6*KS.vz*(no3*KS.vx*KS.vx - no2)*(no3*KS.vy*KS.vy - no2)*KS.fz;
            
        dreal S_0, S_1, S_2, S_3, S_4, S_5, S_6, S_7, S_8, S_9, S_10, S_11, S_12, S_13, S_14, S_15, S_16, S_17, S_18, S_19, S_20, S_21, S_22, S_23, S_24, S_25, S_26;
        
        S_0 = n1o27*(m_0 - no3*m_9 + no3*m_17 - m_26);
        S_1 = n1o108*(no4*m_0 + no6*m_1 + no9*m_7 + no3*m_8 - no6*m_9 - no6*m_10 - no6*m_18 + no6*m_23 + no2*m_26);
        S_2 = n1o108*(no4*m_0 - no6*m_1 + no9*m_7 + no3*m_8 - no6*m_9 + no6*m_10 - no6*m_18 - no6*m_23 + no2*m_26);
        S_3 = n1o108*(no4*m_0 + no6*m_2 - no9*m_7 + no3*m_8 - no6*m_9 - no6*m_11 + no3*m_18 - no9*m_19 + no6*m_24 + no2*m_26);
        S_4 = n1o108*(no4*m_0 - no6*m_2 - no9*m_7 + no3*m_8 - no6*m_9 + no6*m_11 + no3*m_18 - no9*m_19 - no6*m_24 + no2*m_26);
        S_5 = n1o108*(no4*m_0 + no6*m_3 - no6*m_8 - no6*m_9 - no6*m_12 + no3*m_18 + no9*m_19 + no6*m_25 + no2*m_26);
        S_6 = n1o108*(no4*m_0 - no6*m_3 - no6*m_8 - no6*m_9 + no6*m_12 + no3*m_18 + no9*m_19 - no6*m_25 + no2*m_26);
        S_7 = n1o216*(no8*m_0 + no12*m_1 + no12*m_2 + no18*m_4 + no12*m_8 - no3*m_10 - no3*m_11 + no27*m_13 + no27*m_14 - no6*m_17 + no3*m_18 + no9*m_19 - no18*m_22 - no6*m_23 - no6*m_24 - no2*m_26);
        S_8 = n1o216*(no8*m_0 - no12*m_1 + no12*m_2 - no18*m_4 + no12*m_8 + no3*m_10 - no3*m_11 - no27*m_13 + no27*m_14 - no6*m_17 + no3*m_18 + no9*m_19 + no18*m_22 + no6*m_23 - no6*m_24 - no2*m_26);
        S_9 = n1o216*(no8*m_0 + no12*m_1 - no12*m_2 - no18*m_4 + no12*m_8 - no3*m_10 + no3*m_11 + no27*m_13 - no27*m_14 - no6*m_17 + no3*m_18 + no9*m_19 + no18*m_22 - no6*m_23 + no6*m_24 - no2*m_26);
        S_10 = n1o216*(no8*m_0 - no12*m_1 - no12*m_2 + no18*m_4 + no12*m_8 + no3*m_10 + no3*m_11 - no27*m_13 - no27*m_14 - no6*m_17 + no3*m_18 + no9*m_19 - no18*m_22 + no6*m_23 + no6*m_24 - no2*m_26);
        S_11 = n1o216*(no8*m_0 + no12*m_1 + no12*m_3 + no18*m_5 + no18*m_7 - no6*m_8 - no3*m_10 - no3*m_12 - no27*m_13 + no27*m_15 - no6*m_17 + no3*m_18 - no9*m_19 - no18*m_21 - no6*m_23 - no6*m_25 - no2*m_26);
        S_12 = n1o216*(no8*m_0 - no12*m_1 + no12*m_3 - no18*m_5 + no18*m_7 - no6*m_8 + no3*m_10 - no3*m_12 + no27*m_13 + no27*m_15 - no6*m_17 + no3*m_18 - no9*m_19 + no18*m_21 + no6*m_23 - no6*m_25 - no2*m_26);
        S_13 = n1o216*(no8*m_0 + no12*m_1 - no12*m_3 - no18*m_5 + no18*m_7 - no6*m_8 - no3*m_10 + no3*m_12 - no27*m_13 - no27*m_15 - no6*m_17 + no3*m_18 - no9*m_19 + no18*m_21 - no6*m_23 + no6*m_25 - no2*m_26);
	S_14 = n1o216*(no8*m_0 - no12*m_1 - no12*m_3 + no18*m_5 + no18*m_7 - no6*m_8 + no3*m_10 + no3*m_12 + no27*m_13 - no27*m_15 - no6*m_17 + no3*m_18 - no9*m_19 - no18*m_21 + no6*m_23 + no6*m_25 - no2*m_26);
	S_15 = n1o216*(no8*m_0 + no12*m_2 + no12*m_3 + no18*m_6 - no18*m_7 - no6*m_8 - no3*m_11 - no3*m_12 - no27*m_14 - no27*m_15 - no6*m_17 - no6*m_18 - no18*m_20 - no6*m_24 - no6*m_25 - no2*m_26);
        S_16 = n1o216*(no8*m_0 - no12*m_2 + no12*m_3 - no18*m_6 - no18*m_7 - no6*m_8 + no3*m_11 - no3*m_12 + no27*m_14 - no27*m_15 - no6*m_17 - no6*m_18 + no18*m_20 + no6*m_24 - no6*m_25 - no2*m_26);
        S_17 = n1o216*(no8*m_0 + no12*m_2 - no12*m_3 - no18*m_6 - no18*m_7 - no6*m_8 - no3*m_11 + no3*m_12 - no27*m_14 + no27*m_15 - no6*m_17 - no6*m_18 + no18*m_20 - no6*m_24 + no6*m_25 - no2*m_26);
        S_18 = n1o216*(no8*m_0 - no12*m_2 - no12*m_3 + no18*m_6 - no18*m_7 - no6*m_8 + no3*m_11 + no3*m_12 + no27*m_14 + no27*m_15 - no6*m_17 - no6*m_18 - no18*m_20 + no6*m_24 + no6*m_25 - no2*m_26);
        S_19 = n1o216*(no8*m_0 + no12*m_1 + no12*m_2 + no12*m_3 + no18*m_4 + no18*m_5 + no18*m_6 + no12*m_9 + no6*m_10 + no6*m_11 + no6*m_12 + no27*m_16 + no6*m_17 + no9*m_20 + no9*m_21 + no9*m_22 + no3*m_23 + no3*m_24 + no3*m_25 + m_26);
        S_20 = n1o216*(no8*m_0 - no12*m_1 + no12*m_2 + no12*m_3 - no18*m_4 - no18*m_5 + no18*m_6 + no12*m_9 - no6*m_10 + no6*m_11 + no6*m_12 - no27*m_16 + no6*m_17 + no9*m_20 - no9*m_21 - no9*m_22 - no3*m_23 + no3*m_24 + no3*m_25 + m_26);
        S_21 = n1o216*(no8*m_0 + no12*m_1 - no12*m_2 + no12*m_3 - no18*m_4 + no18*m_5 - no18*m_6 + no12*m_9 + no6*m_10 - no6*m_11 + no6*m_12 - no27*m_16 + no6*m_17 - no9*m_20 + no9*m_21 - no9*m_22 + no3*m_23 - no3*m_24 + no3*m_25 + m_26);
        S_22 = n1o216*(no8*m_0 - no12*m_1 - no12*m_2 + no12*m_3 + no18*m_4 - no18*m_5 - no18*m_6 + no12*m_9 - no6*m_10 - no6*m_11 + no6*m_12 + no27*m_16 + no6*m_17 - no9*m_20 - no9*m_21 + no9*m_22 - no3*m_23 - no3*m_24 + no3*m_25 + m_26);
        S_23 = n1o216*(no8*m_0 + no12*m_1 + no12*m_2 - no12*m_3 + no18*m_4 - no18*m_5 - no18*m_6 + no12*m_9 + no6*m_10 + no6*m_11 - no6*m_12 - no27*m_16 + no6*m_17 - no9*m_20 - no9*m_21 + no9*m_22 + no3*m_23 + no3*m_24 - no3*m_25 + m_26);
	S_24 = n1o216*(no8*m_0 - no12*m_1 + no12*m_2 - no12*m_3 - no18*m_4 + no18*m_5 - no18*m_6 + no12*m_9 - no6*m_10 + no6*m_11 - no6*m_12 + no27*m_16 + no6*m_17 - no9*m_20 + no9*m_21 - no9*m_22 - no3*m_23 + no3*m_24 - no3*m_25 + m_26);
	S_25 = n1o216*(no8*m_0 + no12*m_1 - no12*m_2 - no12*m_3 - no18*m_4 - no18*m_5 + no18*m_6 + no12*m_9 + no6*m_10 - no6*m_11 - no6*m_12 + no27*m_16 + no6*m_17 + no9*m_20 - no9*m_21 - no9*m_22 + no3*m_23 - no3*m_24 - no3*m_25 + m_26);
        S_26 = n1o216*(no8*m_0 - no12*m_1 - no12*m_2 - no12*m_3 + no18*m_4 + no18*m_5 + no18*m_6 + no12*m_9 - no6*m_10 - no6*m_11 - no6*m_12 - no27*m_16 + no6*m_17 + no9*m_20 + no9*m_21 + no9*m_22 - no3*m_23 - no3*m_24 - no3*m_25 + m_26);
        
        
        KS.f[zzz] += -(dreal)2.*g_9 + (dreal)4.*g_17 - (dreal)8.*g_26 + S_0; //0
        KS.f[pzz] += g_7 + g_8 - g_9 - (dreal)4.*g_10 - (dreal)4.*g_18 + (dreal)4.*g_23 + (dreal)4.*g_26 + S_1; //1
        KS.f[mzz] += g_7 + g_8 - g_9 + (dreal)4.*g_10 - (dreal)4.*g_18 - (dreal)4.*g_23 + (dreal)4.*g_26 + S_2; //2
        KS.f[zpz] += -g_7 + g_8 - g_9 - (dreal)4.*g_11 + (dreal)2.*g_18 - (dreal)2.*g_19 + (dreal)4.*g_24 + (dreal)4.*g_26 + S_3; //3
        KS.f[zmz] += -g_7 + g_8 - g_9 + (dreal)4.*g_11 + (dreal)2.*g_18 - (dreal)2.*g_19 - (dreal)4.*g_24 + (dreal)4.*g_26 + S_4 ; //4
        KS.f[zzp] += (dreal)-2.*g_8 - g_9 - (dreal)4.*g_12 + (dreal)2.*g_18 + (dreal)2.*g_19 + (dreal)4.*g_25 + (dreal)4.*g_26 + S_5; //5
        KS.f[zzm] += (dreal)-2.*g_8 - g_9 + (dreal)4.*g_12 + (dreal)2.*g_18 + (dreal)2.*g_19 - (dreal)4.*g_25 + (dreal)4.*g_26 + S_6; //6
        KS.f[ppz] += g_4 + (dreal)2.*g_8 - g_10 - g_11 + g_13 + g_14 - g_17 + g_18 + g_19 - (dreal)2.*g_22 - (dreal)2.*g_23 - (dreal)2.*g_24 - (dreal)2.*g_26 + S_7; //7
        KS.f[mpz] += -g_4 + (dreal)2.*g_8 + g_10 - g_11 - g_13 + g_14 - g_17 + g_18 + g_19 + (dreal)2.*g_22 + (dreal)2.*g_23 - (dreal)2.*g_24 - (dreal)2.*g_26 + S_8; //8
        KS.f[pmz] += -g_4 + (dreal)2.*g_8 - g_10 + g_11 + g_13 - g_14 - g_17 + g_18 + g_19 + (dreal)2.*g_22 - (dreal)2.*g_23 + (dreal)2.*g_24 - (dreal)2.*g_26 + S_9; //9
        KS.f[mmz] += g_4 + (dreal)2.*g_8 + g_10 + g_11 - g_13 - g_14 - g_17 + g_18 + g_19 - (dreal)2.*g_22 + (dreal)2.*g_23 + (dreal)2.*g_24 - (dreal)2.*g_26 + S_10; //10
        KS.f[pzp] += g_5 + g_7 - g_8 - g_10 - g_12 - g_13 + g_15 - g_17 + g_18 - g_19 - (dreal)2.*g_21 - (dreal)2.*g_23 - (dreal)2.*g_25 - (dreal)2.*g_26 + S_11; //11
        KS.f[mzp] += -g_5 + g_7 - g_8 + g_10 - g_12 + g_13 + g_15 - g_17 + g_18 - g_19 + (dreal)2.*g_21 + (dreal)2.*g_23 - (dreal)2.*g_25 - (dreal)2.*g_26 + S_12; //12
        KS.f[pzm] += -g_5 + g_7 - g_8 - g_10 + g_12 - g_13 - g_15 - g_17 + g_18 - g_19 + (dreal)2.*g_21 - (dreal)2.*g_23 + (dreal)2.*g_25 - (dreal)2.*g_26 + S_13; //13
        KS.f[mzm] += g_5 + g_7 - g_8 + g_10 + g_12 + g_13 - g_15 - g_17 + g_18 - g_19 - (dreal)2.*g_21 + (dreal)2.*g_23 + (dreal)2.*g_25 - (dreal)2.*g_26 + S_14; //14
        KS.f[zpp] += g_6 - g_7 - g_8 - g_11 - g_12 - g_14 - g_15 - g_17 - (dreal)2.*g_18 - (dreal)2.*g_20 - (dreal)2.*g_24 - (dreal)2.*g_25 - (dreal)2.*g_26 + S_15; //15
        KS.f[zmp] += -g_6 - g_7 - g_8 + g_11 - g_12 + g_14 - g_15 - g_17 - (dreal)2.*g_18 + (dreal)2.*g_20 + (dreal)2.*g_24 - (dreal)2.*g_25 - (dreal)2.*g_26 + S_16; //16
        KS.f[zpm] += -g_6 - g_7 - g_8 - g_11 + g_12 - g_14 + g_15 - g_17 - (dreal)2.*g_18 + (dreal)2.*g_20 - (dreal)2.*g_24 + (dreal)2.*g_25 - (dreal)2.*g_26 + S_17; //17
        KS.f[zmm] += g_6 - g_7 - g_8 + g_11 + g_12 + g_14 + g_15 - g_17 - (dreal)2.*g_18 - (dreal)2.*g_20 + (dreal)2.*g_24 + (dreal)2.*g_25 - (dreal)2.*g_26 + S_18; //18
        KS.f[ppp] += g_4 + g_5 + g_6 + g_9 + (dreal)2.*g_10 + (dreal)2.*g_11 + (dreal)2.*g_12 + g_16 + g_17 + g_20 + g_21 + g_22 + g_23 + g_24 + g_25 + g_26 + S_19; //19
        KS.f[mpp] += -g_4 - g_5 + g_6 + g_9 - (dreal)2.*g_10 + (dreal)2.*g_11 + (dreal)2.*g_12 - g_16 + g_17 + g_20 - g_21 - g_22 - g_23 + g_24 + g_25 + g_26 + S_20; //20
        KS.f[pmp] += -g_4 + g_5 - g_6 + g_9 + (dreal)2.*g_10 - (dreal)2.*g_11 + (dreal)2.*g_12 - g_16 + g_17 - g_20 + g_21 - g_22 + g_23 - g_24 + g_25 + g_26 + S_21; //21
        KS.f[mmp] += g_4 - g_5 - g_6 + g_9 - (dreal)2.*g_10 - (dreal)2.*g_11 + (dreal)2.*g_12 + g_16 + g_17 - g_20 - g_21 + g_22 - g_23 - g_24 + g_25 + g_26 + S_22; //22
        KS.f[ppm] += g_4 - g_5 - g_6 + g_9 + (dreal)2.*g_10 + (dreal)2.*g_11 - (dreal)2.*g_12 - g_16 + g_17 - g_20 - g_21 + g_22 + g_23 + g_24 - g_25 + g_26 + S_23; //23
        KS.f[mpm] += -g_4 + g_5 - g_6 + g_9 - (dreal)2.*g_10 + (dreal)2.*g_11 - (dreal)2.*g_12 + g_16 + g_17 - g_20 + g_21 - g_22 - g_23 + g_24 - g_25 + g_26 + S_24; //24
        KS.f[pmm] += -g_4 - g_5 + g_6 + g_9 + (dreal)2.*g_10 - (dreal)2.*g_11 - (dreal)2.*g_12 + g_16 + g_17 + g_20 - g_21 - g_22 + g_23 - g_24 - g_25 + g_26 + S_25; //25
        KS.f[mmm] += g_4 + g_5 + g_6 + g_9 - (dreal)2.*g_10 - (dreal)2.*g_11 - (dreal)2.*g_12 - g_16 + g_17 + g_20 + g_21 + g_22 - g_23 - g_24 - g_25 + g_26 + S_26; //26
         
        
	}
};

