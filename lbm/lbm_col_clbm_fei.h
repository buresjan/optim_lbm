#include "lbm_common.h"
// FEI and LUO standard CLBM ... not working
template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_CLBM_FEI : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "CLBM_FEI"; 
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		// CLBM Fei & Luo & Li
		// generated from Maple file: clbm3d_fei_luo_li.mw
		// step 1: raw moments 
		const dreal t1 = KS.rho; //KS.f[zzz]+KS.f[pzz]+KS.f[mzz]+KS.f[zpz]+KS.f[zmz]+KS.f[zzp]+KS.f[zzm]+KS.f[ppz]+KS.f[mpz]+KS.f[pmz]+KS.f[mmz]+KS.f[pzp]+KS.f[mzp]+KS.f[pzm]+KS.f[mzm]+KS.f[zpp]+KS.f[zmp]+KS.f[zpm]+KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm];
		const dreal t2 = KS.f[pzz]-KS.f[mzz]+KS.f[ppz]-KS.f[mpz]+KS.f[pmz]-KS.f[mmz]+KS.f[pzp]-KS.f[mzp]+KS.f[pzm]-KS.f[mzm]+KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]+KS.f[ppm]-KS.f[mpm]+KS.f[pmm]-KS.f[mmm];
		const dreal t3 = KS.f[zpz]-KS.f[zmz]+KS.f[ppz]+KS.f[mpz]-KS.f[pmz]-KS.f[mmz]+KS.f[zpp]-KS.f[zmp]+KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]+KS.f[ppm]+KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
		const dreal t4 = KS.f[zzp]-KS.f[zzm]+KS.f[pzp]+KS.f[mzp]-KS.f[pzm]-KS.f[mzm]+KS.f[zpp]+KS.f[zmp]-KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]-KS.f[ppm]-KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
//		const dreal t2 = KS.rho*KS.vx; //KS.f[pzz]-KS.f[mzz]+KS.f[ppz]-KS.f[mpz]+KS.f[pmz]-KS.f[mmz]+KS.f[pzp]-KS.f[mzp]+KS.f[pzm]-KS.f[mzm]+KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]+KS.f[ppm]-KS.f[mpm]+KS.f[pmm]-KS.f[mmm];
//		const dreal t3 = KS.rho*KS.vy; //KS.f[zpz]-KS.f[zmz]+KS.f[ppz]+KS.f[mpz]-KS.f[pmz]-KS.f[mmz]+KS.f[zpp]-KS.f[zmp]+KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]+KS.f[ppm]+KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
//		const dreal t4 = KS.rho*KS.vz; //KS.f[zzp]-KS.f[zzm]+KS.f[pzp]+KS.f[mzp]-KS.f[pzm]-KS.f[mzm]+KS.f[zpp]+KS.f[zmp]-KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]-KS.f[ppm]-KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
		const dreal t5 = KS.f[ppz]-KS.f[mpz]-KS.f[pmz]+KS.f[mmz]+KS.f[ppp]-KS.f[mpp]-KS.f[pmp]+KS.f[mmp]+KS.f[ppm]-KS.f[mpm]-KS.f[pmm]+KS.f[mmm];
		const dreal t6 = KS.f[pzp]-KS.f[mzp]-KS.f[pzm]+KS.f[mzm]+KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]-KS.f[ppm]+KS.f[mpm]-KS.f[pmm]+KS.f[mmm];
		const dreal t7 = KS.f[zpp]-KS.f[zmp]-KS.f[zpm]+KS.f[zmm]+KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]-KS.f[ppm]-KS.f[mpm]+KS.f[pmm]+KS.f[mmm];
		const dreal t8 = KS.f[pzz]+KS.f[mzz]+KS.f[ppz]+KS.f[mpz]+KS.f[pmz]+KS.f[mmz]+KS.f[pzp]+KS.f[mzp]+KS.f[pzm]+KS.f[mzm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm];
		const dreal t9 = KS.f[zpz]+KS.f[zmz]+KS.f[ppz]+KS.f[mpz]+KS.f[pmz]+KS.f[mmz]+KS.f[zpp]+KS.f[zmp]+KS.f[zpm]+KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm];
		const dreal t10 = KS.f[zzp]+KS.f[zzm]+KS.f[pzp]+KS.f[mzp]+KS.f[pzm]+KS.f[mzm]+KS.f[zpp]+KS.f[zmp]+KS.f[zpm]+KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm]; 
		const dreal t11 = KS.f[ppz]-KS.f[mpz]+KS.f[pmz]-KS.f[mmz]+KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]+KS.f[ppm]-KS.f[mpm]+KS.f[pmm]-KS.f[mmm]; 
		const dreal t12 = KS.f[pzp]-KS.f[mzp]+KS.f[pzm]-KS.f[mzm]+KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]+KS.f[ppm]-KS.f[mpm]+KS.f[pmm]-KS.f[mmm];
		const dreal t13 = KS.f[ppz]+KS.f[mpz]-KS.f[pmz]-KS.f[mmz]+KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]+KS.f[ppm]+KS.f[mpm]-KS.f[pmm]-KS.f[mmm]; 
		const dreal t14 = KS.f[pzp]+KS.f[mzp]-KS.f[pzm]-KS.f[mzm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]-KS.f[ppm]-KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
		const dreal t15 = KS.f[zpp]-KS.f[zmp]+KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]+KS.f[ppm]+KS.f[mpm]-KS.f[pmm]-KS.f[mmm]; 
		const dreal t16 = KS.f[zpp]+KS.f[zmp]-KS.f[zpm]-KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]-KS.f[ppm]-KS.f[mpm]-KS.f[pmm]-KS.f[mmm]; 
		const dreal t17 = KS.f[ppp]-KS.f[mpp]-KS.f[pmp]+KS.f[mmp]-KS.f[ppm]+KS.f[mpm]+KS.f[pmm]-KS.f[mmm]; 
		const dreal t18 = KS.f[ppz]+KS.f[mpz]+KS.f[pmz]+KS.f[mmz]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm];
		const dreal t19 = KS.f[pzp]+KS.f[mzp]+KS.f[pzm]+KS.f[mzm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm]; 
		const dreal t20 = KS.f[zpp]+KS.f[zmp]+KS.f[zpm]+KS.f[zmm]+KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm]; 
		const dreal t21 = KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]-KS.f[ppm]-KS.f[mpm]+KS.f[pmm]+KS.f[mmm]; 
		const dreal t22 = KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]-KS.f[ppm]+KS.f[mpm]-KS.f[pmm]+KS.f[mmm]; 
		const dreal t23 = KS.f[ppp]-KS.f[mpp]-KS.f[pmp]+KS.f[mmp]+KS.f[ppm]-KS.f[mpm]-KS.f[pmm]+KS.f[mmm]; 
		const dreal t24 = KS.f[ppp]-KS.f[mpp]+KS.f[pmp]-KS.f[mmp]+KS.f[ppm]-KS.f[mpm]+KS.f[pmm]-KS.f[mmm]; 
		const dreal t25 = KS.f[ppp]+KS.f[mpp]-KS.f[pmp]-KS.f[mmp]+KS.f[ppm]+KS.f[mpm]-KS.f[pmm]-KS.f[mmm];
		const dreal t26 = KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]-KS.f[ppm]-KS.f[mpm]-KS.f[pmm]-KS.f[mmm]; 
		const dreal t27 = KS.f[ppp]+KS.f[mpp]+KS.f[pmp]+KS.f[mmp]+KS.f[ppm]+KS.f[mpm]+KS.f[pmm]+KS.f[mmm];

		// step 3: apply S (relaxation)
		const dreal s_1 = no1;///(no3*KS.lbmViscosity+n1o2);
		const dreal s_nu = no1;///(no3*KS.lbmViscosity+n1o2);
//		const dreal s_2 = no1;//no1;//no1/(no3*KS.lbmViscosity+n1o2);
//		const dreal s_2b = s_nu; // FIXME s_2b = no1/(n9o2 * xi + n1o2);
		const dreal s_m = 0;//n1o3*(s_2b-s_2);
		const dreal s_p = s_nu;//n1o3*(s_2b+no2*s_2);
		const dreal s_3 = no1;//(no16-no8*s_2)/(no8-s_2);
		const dreal s_3b= s_3;
		const dreal s_4 = no1;
		const dreal s_4b = no1;
		const dreal s_5 = no1;
		const dreal s_6 = no1;

		// a = tildeT + C + S [ -tildeT + tildeT^eq - 1/2 C ]
		const dreal a1 = t1;
		const dreal a2 = -(n1o2)*s_1*KS.fx+KS.fx+t2;
		const dreal a3 = -(n1o2)*s_1*KS.fy+KS.fy+t3;
		const dreal a4 = -(n1o2)*s_1*KS.fz+KS.fz+t4;
		const dreal a5 = s_nu*(KS.rho*KS.vx*KS.vy-t5)+t5;
		const dreal a6 = s_nu*(KS.rho*KS.vx*KS.vz-t6)+t6;
		const dreal a7 = s_nu*(KS.rho*KS.vy*KS.vz-t7)+t7;
		const dreal a8 = s_m*((n1o3)*KS.rho+KS.vy*KS.vy*KS.rho-t9)+s_m*((n1o3)*KS.rho+KS.vz*KS.vz*KS.rho-t10)+s_p*((n1o3)*KS.rho+KS.vx*KS.vx*KS.rho-t8)+t8;
		const dreal a9 = s_m*((n1o3)*KS.rho+KS.vx*KS.vx*KS.rho-t8)+s_m*((n1o3)*KS.rho+KS.vz*KS.vz*KS.rho-t10)+s_p*((n1o3)*KS.rho+KS.vy*KS.vy*KS.rho-t9)+t9;
		const dreal a10 = s_m*((n1o3)*KS.rho+KS.vx*KS.vx*KS.rho-t8)+s_m*((n1o3)*KS.rho+KS.vy*KS.vy*KS.rho-t9)+s_p*((n1o3)*KS.rho+KS.vz*KS.vz*KS.rho-t10)+t10;
		const dreal a11 = t11+(n1o3)*KS.fx+s_3*(-no2*KS.vx*KS.vy*KS.vy*KS.rho+no2*KS.vy*t5+KS.vx*t9-t11-(n1o6)*KS.fx);
		const dreal a12 = t12+(n1o3)*KS.fx+s_3*(-no2*KS.vx*KS.vz*KS.vz*KS.rho+no2*KS.vz*t6+KS.vx*t10-t12-(n1o6)*KS.fx);
		const dreal a13 = t13+(n1o3)*KS.fy+s_3*(-no2*KS.vx*KS.vx*KS.vy*KS.rho+no2*KS.vx*t5+KS.vy*t8-t13-(n1o6)*KS.fy);
		const dreal a14 = t14+(n1o3)*KS.fz+s_3*(-no2*KS.vx*KS.vx*KS.vz*KS.rho+no2*KS.vx*t6+KS.vz*t8-t14-(n1o6)*KS.fz);
		const dreal a15 = t15+(n1o3)*KS.fy+s_3*(-no2*KS.vy*KS.vz*KS.vz*KS.rho+no2*KS.vz*t7+KS.vy*t10-t15-(n1o6)*KS.fy);
		const dreal a16 = t16+(n1o3)*KS.fz+s_3*(-no2*KS.vy*KS.vy*KS.vz*KS.rho+no2*KS.vy*t7+KS.vz*t9-t16-(n1o6)*KS.fz);
		const dreal a17 = s_3b*(-no2*KS.rho*KS.vx*KS.vy*KS.vz+t5*KS.vz+t6*KS.vy+t7*KS.vx-t17)+t17;
		const dreal a18 = s_4*((n1o9)*KS.rho+no3*KS.vx*KS.vx*KS.vy*KS.vy*KS.rho-no4*KS.vx*KS.vy*t5-KS.vy*KS.vy*t8-KS.vx*KS.vx*t9+no2*KS.vx*t11+no2*KS.vy*t13-t18)+t18;
		const dreal a19 = s_4*((n1o9)*KS.rho+no3*KS.vx*KS.vx*KS.vz*KS.vz*KS.rho-no4*KS.vx*KS.vz*t6-KS.vz*KS.vz*t8-KS.vx*KS.vx*t10+no2*KS.vx*t12+no2*KS.vz*t14-t19)+t19;
		const dreal a20 = s_4*((n1o9)*KS.rho+no3*KS.vy*KS.vy*KS.vz*KS.vz*KS.rho-no4*KS.vy*KS.vz*t7-KS.vz*KS.vz*t9-KS.vy*KS.vy*t10+no2*KS.vy*t15+no2*KS.vz*t16-t20)+t20;
		const dreal a21 = s_4b*(no3*KS.rho*KS.vx*KS.vx*KS.vy*KS.vz-no2*t5*KS.vx*KS.vz-no2*t6*KS.vx*KS.vy-t7*KS.vx*KS.vx-t8*KS.vy*KS.vz+t13*KS.vz+t14*KS.vy+no2*t17*KS.vx-t21)+t21;
		const dreal a22 = s_4b*(no3*KS.rho*KS.vx*KS.vy*KS.vy*KS.vz-no2*t5*KS.vy*KS.vz-t6*KS.vy*KS.vy-no2*t7*KS.vx*KS.vy-t9*KS.vx*KS.vz+t11*KS.vz+t16*KS.vx+no2*t17*KS.vy-t22)+t22;
		const dreal a23 = s_4b*(no3*KS.rho*KS.vx*KS.vy*KS.vz*KS.vz-t5*KS.vz*KS.vz-no2*t6*KS.vy*KS.vz-no2*t7*KS.vx*KS.vz-t10*KS.vx*KS.vy+t12*KS.vy+t15*KS.vx+no2*t17*KS.vz-t23)+t23;
		const dreal a24 = t24+(n1o9)*KS.fx+s_5*(-no4*KS.vx*KS.vy*KS.vy*KS.vz*KS.vz*KS.rho+no2*KS.vy*KS.vz*KS.vz*t5+no2*KS.vy*KS.vy*KS.vz*t6+no4*KS.vx*KS.vy*KS.vz*t7+KS.vx*KS.vz*KS.vz*t9+KS.vx*KS.vy*KS.vy*t10-KS.vz*KS.vz*t11-KS.vy*KS.vy*t12-no2*KS.vx*KS.vy*t15-no2*KS.vx*KS.vz*t16-no4*KS.vy*KS.vz*t17+KS.vx*t20+no2*KS.vz*t22+no2*KS.vy*t23-t24-(n1o18)*KS.fx);
		const dreal a25 = t25+(n1o9)*KS.fy+s_5*(-no4*KS.vx*KS.vx*KS.vy*KS.vz*KS.vz*KS.rho+no2*KS.vx*KS.vz*KS.vz*t5+no4*KS.vx*KS.vy*KS.vz*t6+no2*KS.vx*KS.vx*KS.vz*t7+KS.vy*KS.vz*KS.vz*t8+KS.vx*KS.vx*KS.vy*t10-no2*KS.vx*KS.vy*t12-KS.vz*KS.vz*t13-no2*KS.vy*KS.vz*t14-KS.vx*KS.vx*t15-no4*KS.vx*KS.vz*t17+KS.vy*t19+no2*KS.vz*t21+no2*KS.vx*t23-t25-(n1o18)*KS.fy);
		const dreal a26 = t26+(n1o9)*KS.fz+s_5*(-no4*KS.vx*KS.vx*KS.vy*KS.vy*KS.vz*KS.rho+no4*KS.vx*KS.vy*KS.vz*t5+no2*KS.vx*KS.vy*KS.vy*t6+no2*KS.vx*KS.vx*KS.vy*t7+KS.vy*KS.vy*KS.vz*t8+KS.vx*KS.vx*KS.vz*t9-no2*KS.vx*KS.vz*t11-no2*KS.vy*KS.vz*t13-KS.vy*KS.vy*t14-KS.vx*KS.vx*t16-no4*KS.vx*KS.vy*t17+KS.vz*t18+no2*KS.vy*t21+no2*KS.vx*t22-t26-(n1o18)*KS.fz);
		const dreal a27 = s_6*(-KS.vy*KS.vy*KS.vz*KS.vz*t8-KS.vx*KS.vx*KS.vz*KS.vz*t9-KS.vx*KS.vx*KS.vy*KS.vy*t10+no2*KS.vx*KS.vz*KS.vz*t11+no2*KS.vx*KS.vy*KS.vy*t12+no2*KS.vy*KS.vz*KS.vz*t13+no2*KS.vy*KS.vy*KS.vz*t14+no2*KS.vx*KS.vx*KS.vy*t15+no2*KS.vx*KS.vx*KS.vz*t16-no4*KS.vy*KS.vz*t21-no4*KS.vx*KS.vz*t22-no4*KS.vx*KS.vy*t23-t27+(n1o27)*KS.rho-KS.vz*KS.vz*t18-KS.vy*KS.vy*t19-KS.vx*KS.vx*t20+no2*KS.vx*t24+no2*KS.vy*t25+no2*KS.vz*t26-no4*KS.vx*KS.vy*KS.vz*KS.vz*t5-no4*KS.vx*KS.vy*KS.vy*KS.vz*t6-no4*KS.vx*KS.vx*KS.vy*KS.vz*t7+no8*KS.vx*KS.vy*KS.vz*t17+no5*KS.vx*KS.vx*KS.vy*KS.vy*KS.vz*KS.vz*KS.rho)+t27;

		const dreal b1 = a1;
		const dreal b2 = a1*KS.vx+a2;
		const dreal b3 = a1*KS.vy+a3;
		const dreal b4 = a1*KS.vz+a4;
		const dreal b5 = a1*KS.vx*KS.vy+a2*KS.vy+a3*KS.vx+a5;
		const dreal b6 = a1*KS.vx*KS.vz+a2*KS.vz+a4*KS.vx+a6;
		const dreal b7 = a1*KS.vy*KS.vz+a3*KS.vz+a4*KS.vy+a7;
		const dreal b8 = a1*KS.vx*KS.vx+no2*a2*KS.vx+a8;
		const dreal b9 = a1*KS.vy*KS.vy+no2*a3*KS.vy+a9;
		const dreal b10 = a1*KS.vz*KS.vz+no2*a4*KS.vz+a10;
		const dreal b11 = a1*KS.vx*KS.vy*KS.vy+a2*KS.vy*KS.vy+no2*a3*KS.vx*KS.vy+no2*a5*KS.vy+a9*KS.vx+a11;
		const dreal b12 = a1*KS.vx*KS.vz*KS.vz+a2*KS.vz*KS.vz+no2*a4*KS.vx*KS.vz+no2*a6*KS.vz+a10*KS.vx+a12;
		const dreal b13 = a1*KS.vx*KS.vx*KS.vy+no2*a2*KS.vx*KS.vy+a3*KS.vx*KS.vx+no2*a5*KS.vx+a8*KS.vy+a13;
		const dreal b14 = a1*KS.vx*KS.vx*KS.vz+no2*a2*KS.vx*KS.vz+a4*KS.vx*KS.vx+no2*a6*KS.vx+a8*KS.vz+a14;
		const dreal b15 = a1*KS.vy*KS.vz*KS.vz+a3*KS.vz*KS.vz+no2*a4*KS.vy*KS.vz+no2*a7*KS.vz+a10*KS.vy+a15;
		const dreal b16 = a1*KS.vy*KS.vy*KS.vz+no2*a3*KS.vy*KS.vz+a4*KS.vy*KS.vy+no2*a7*KS.vy+a9*KS.vz+a16;
		const dreal b17 = a1*KS.vx*KS.vy*KS.vz+a2*KS.vy*KS.vz+a3*KS.vx*KS.vz+a4*KS.vx*KS.vy+a5*KS.vz+a6*KS.vy+a7*KS.vx+a17;
		const dreal b18 = a1*KS.vx*KS.vx*KS.vy*KS.vy+no2*a2*KS.vx*KS.vy*KS.vy+no2*a3*KS.vx*KS.vx*KS.vy+no4*a5*KS.vx*KS.vy+a8*KS.vy*KS.vy+a9*KS.vx*KS.vx+no2*a11*KS.vx+no2*a13*KS.vy+a18;
		const dreal b19 = a1*KS.vx*KS.vx*KS.vz*KS.vz+no2*a2*KS.vx*KS.vz*KS.vz+no2*a4*KS.vx*KS.vx*KS.vz+no4*a6*KS.vx*KS.vz+a8*KS.vz*KS.vz+a10*KS.vx*KS.vx+no2*a12*KS.vx+no2*a14*KS.vz+a19;
		const dreal b20 = a1*KS.vy*KS.vy*KS.vz*KS.vz+no2*a3*KS.vy*KS.vz*KS.vz+no2*a4*KS.vy*KS.vy*KS.vz+no4*a7*KS.vy*KS.vz+a9*KS.vz*KS.vz+a10*KS.vy*KS.vy+no2*a15*KS.vy+no2*a16*KS.vz+a20;
		const dreal b21 = a1*KS.vx*KS.vx*KS.vy*KS.vz+no2*a2*KS.vx*KS.vy*KS.vz+a3*KS.vx*KS.vx*KS.vz+a4*KS.vx*KS.vx*KS.vy+no2*a5*KS.vx*KS.vz+no2*a6*KS.vx*KS.vy+a7*KS.vx*KS.vx+a8*KS.vy*KS.vz+a13*KS.vz+a14*KS.vy+no2*a17*KS.vx+a21;
		const dreal b22 = a1*KS.vx*KS.vy*KS.vy*KS.vz+a2*KS.vy*KS.vy*KS.vz+no2*a3*KS.vx*KS.vy*KS.vz+a4*KS.vx*KS.vy*KS.vy+no2*a5*KS.vy*KS.vz+a6*KS.vy*KS.vy+no2*a7*KS.vx*KS.vy+a9*KS.vx*KS.vz+a11*KS.vz+a16*KS.vx+no2*a17*KS.vy+a22;
		const dreal b23 = a1*KS.vx*KS.vy*KS.vz*KS.vz+a2*KS.vy*KS.vz*KS.vz+a3*KS.vx*KS.vz*KS.vz+no2*a4*KS.vx*KS.vy*KS.vz+a5*KS.vz*KS.vz+no2*a6*KS.vy*KS.vz+no2*a7*KS.vx*KS.vz+a10*KS.vx*KS.vy+a12*KS.vy+a15*KS.vx+no2*a17*KS.vz+a23;
		const dreal b24 = a1*KS.vx*KS.vy*KS.vy*KS.vz*KS.vz+a2*KS.vy*KS.vy*KS.vz*KS.vz+no2*a3*KS.vx*KS.vy*KS.vz*KS.vz+no2*a4*KS.vx*KS.vy*KS.vy*KS.vz+no2*a5*KS.vy*KS.vz*KS.vz+no2*a6*KS.vy*KS.vy*KS.vz+no4*a7*KS.vx*KS.vy*KS.vz+a9*KS.vx*KS.vz*KS.vz+a10*KS.vx*KS.vy*KS.vy+a11*KS.vz*KS.vz+a12*KS.vy*KS.vy+no2*a15*KS.vx*KS.vy+no2*a16*KS.vx*KS.vz+no4*a17*KS.vy*KS.vz+a20*KS.vx+no2*a22*KS.vz+no2*a23*KS.vy+a24;
		const dreal b25 = a1*KS.vx*KS.vx*KS.vy*KS.vz*KS.vz+no2*a2*KS.vx*KS.vy*KS.vz*KS.vz+a3*KS.vx*KS.vx*KS.vz*KS.vz+no2*a4*KS.vx*KS.vx*KS.vy*KS.vz+no2*a5*KS.vx*KS.vz*KS.vz+no4*a6*KS.vx*KS.vy*KS.vz+no2*a7*KS.vx*KS.vx*KS.vz+a8*KS.vy*KS.vz*KS.vz+a10*KS.vx*KS.vx*KS.vy+no2*a12*KS.vx*KS.vy+a13*KS.vz*KS.vz+no2*a14*KS.vy*KS.vz+a15*KS.vx*KS.vx+no4*a17*KS.vx*KS.vz+a19*KS.vy+no2*a21*KS.vz+no2*a23*KS.vx+a25;
		const dreal b26 = a1*KS.vx*KS.vx*KS.vy*KS.vy*KS.vz+no2*a2*KS.vx*KS.vy*KS.vy*KS.vz+no2*a3*KS.vx*KS.vx*KS.vy*KS.vz+a4*KS.vx*KS.vx*KS.vy*KS.vy+no4*a5*KS.vx*KS.vy*KS.vz+no2*a6*KS.vx*KS.vy*KS.vy+no2*a7*KS.vx*KS.vx*KS.vy+a8*KS.vy*KS.vy*KS.vz+a9*KS.vx*KS.vx*KS.vz+no2*a11*KS.vx*KS.vz+no2*a13*KS.vy*KS.vz+a14*KS.vy*KS.vy+a16*KS.vx*KS.vx+no4*a17*KS.vx*KS.vy+a18*KS.vz+no2*a21*KS.vy+no2*a22*KS.vx+a26;
		const dreal b27 = a1*KS.vx*KS.vx*KS.vy*KS.vy*KS.vz*KS.vz+no2*a2*KS.vx*KS.vy*KS.vy*KS.vz*KS.vz+no2*a3*KS.vx*KS.vx*KS.vy*KS.vz*KS.vz+no2*a4*KS.vx*KS.vx*KS.vy*KS.vy*KS.vz+no4*a5*KS.vx*KS.vy*KS.vz*KS.vz+no4*a6*KS.vx*KS.vy*KS.vy*KS.vz+no4*a7*KS.vx*KS.vx*KS.vy*KS.vz+a8*KS.vy*KS.vy*KS.vz*KS.vz+a9*KS.vx*KS.vx*KS.vz*KS.vz+a10*KS.vx*KS.vx*KS.vy*KS.vy+no2*a11*KS.vx*KS.vz*KS.vz+no2*a12*KS.vx*KS.vy*KS.vy+no2*a13*KS.vy*KS.vz*KS.vz+no2*a14*KS.vy*KS.vy*KS.vz+no2*a15*KS.vx*KS.vx*KS.vy+no2*a16*KS.vx*KS.vx*KS.vz+no8*a17*KS.vx*KS.vy*KS.vz+a18*KS.vz*KS.vz+a19*KS.vy*KS.vy+a20*KS.vx*KS.vx+no4*a21*KS.vy*KS.vz+no4*a22*KS.vx*KS.vz+no4*a23*KS.vx*KS.vy+no2*a24*KS.vx+no2*a25*KS.vy+no2*a26*KS.vz+a27;

		KS.f[zzz] = b1-b8-b9-b10+b18+b19+b20-b27;
		KS.f[pzz] = n1o2*b2+n1o2*b8-n1o2*b11-n1o2*b12-n1o2*b18-n1o2*b19+n1o2*b24+n1o2*b27;
		KS.f[mzz] = -n1o2*b2+n1o2*b8+n1o2*b11+n1o2*b12-n1o2*b18-n1o2*b19-n1o2*b24+n1o2*b27;
		KS.f[zpz] = n1o2*b3+n1o2*b9-n1o2*b13-n1o2*b15-n1o2*b18-n1o2*b20+n1o2*b25+n1o2*b27;
		KS.f[zmz] = -n1o2*b3+n1o2*b9+n1o2*b13+n1o2*b15-n1o2*b18-n1o2*b20-n1o2*b25+n1o2*b27;
		KS.f[zzp] = n1o2*b4+n1o2*b10-n1o2*b14-n1o2*b16-n1o2*b19-n1o2*b20+n1o2*b26+n1o2*b27;
		KS.f[zzm] = -n1o2*b4+n1o2*b10+n1o2*b14+n1o2*b16-n1o2*b19-n1o2*b20-n1o2*b26+n1o2*b27;
		KS.f[ppz] = n1o4*b5+n1o4*b11+n1o4*b13+n1o4*b18-n1o4*b23-n1o4*b24-n1o4*b25-n1o4*b27;
		KS.f[mpz] = -n1o4*b5-n1o4*b11+n1o4*b13+n1o4*b18+n1o4*b23+n1o4*b24-n1o4*b25-n1o4*b27;
		KS.f[pmz] = -n1o4*b5+n1o4*b11-n1o4*b13+n1o4*b18+n1o4*b23-n1o4*b24+n1o4*b25-n1o4*b27;
		KS.f[mmz] = n1o4*b5-n1o4*b11-n1o4*b13+n1o4*b18-n1o4*b23+n1o4*b24+n1o4*b25-n1o4*b27;
		KS.f[pzp] = n1o4*b6+n1o4*b12+n1o4*b14+n1o4*b19-n1o4*b22-n1o4*b24-n1o4*b26-n1o4*b27;
		KS.f[mzp] = -n1o4*b6-n1o4*b12+n1o4*b14+n1o4*b19+n1o4*b22+n1o4*b24-n1o4*b26-n1o4*b27;
		KS.f[pzm] = -n1o4*b6+n1o4*b12-n1o4*b14+n1o4*b19+n1o4*b22-n1o4*b24+n1o4*b26-n1o4*b27;
		KS.f[mzm] = n1o4*b6-n1o4*b12-n1o4*b14+n1o4*b19-n1o4*b22+n1o4*b24+n1o4*b26-n1o4*b27;
		KS.f[zpp] = n1o4*b7+n1o4*b15+n1o4*b16+n1o4*b20-n1o4*b21-n1o4*b25-n1o4*b26-n1o4*b27;
		KS.f[zmp] = -n1o4*b7-n1o4*b15+n1o4*b16+n1o4*b20+n1o4*b21+n1o4*b25-n1o4*b26-n1o4*b27;
		KS.f[zpm] = -n1o4*b7+n1o4*b15-n1o4*b16+n1o4*b20+n1o4*b21-n1o4*b25+n1o4*b26-n1o4*b27;
		KS.f[zmm] = n1o4*b7-n1o4*b15-n1o4*b16+n1o4*b20-n1o4*b21+n1o4*b25+n1o4*b26-n1o4*b27;
		KS.f[ppp] = n1o8*b17+n1o8*b21+n1o8*b22+n1o8*b23+n1o8*b24+n1o8*b25+n1o8*b26+n1o8*b27;
		KS.f[mpp] = -n1o8*b17+n1o8*b21-n1o8*b22-n1o8*b23-n1o8*b24+n1o8*b25+n1o8*b26+n1o8*b27;
		KS.f[pmp] = -n1o8*b17-n1o8*b21+n1o8*b22-n1o8*b23+n1o8*b24-n1o8*b25+n1o8*b26+n1o8*b27;
		KS.f[mmp] = n1o8*b17-n1o8*b21-n1o8*b22+n1o8*b23-n1o8*b24-n1o8*b25+n1o8*b26+n1o8*b27;
		KS.f[ppm] = -n1o8*b17-n1o8*b21-n1o8*b22+n1o8*b23+n1o8*b24+n1o8*b25-n1o8*b26+n1o8*b27;
		KS.f[mpm] = n1o8*b17-n1o8*b21+n1o8*b22-n1o8*b23-n1o8*b24+n1o8*b25-n1o8*b26+n1o8*b27;
		KS.f[pmm] = n1o8*b17+n1o8*b21-n1o8*b22-n1o8*b23+n1o8*b24-n1o8*b25-n1o8*b26+n1o8*b27;
		KS.f[mmm] = -n1o8*b17+n1o8*b21+n1o8*b22+n1o8*b23-n1o8*b24-n1o8*b25-n1o8*b26+n1o8*b27;
	}
};
